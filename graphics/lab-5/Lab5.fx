//--------------------------------------------------------------------------------------
// File: Lab5.fx
//
// Lab 5. Author: Andreev Daniil
// Based on: Tutorial07.fx
//
// The textured cube is lit by a point light source:
// ambient + diffuse + specular with distance attenuation.
//
// Copyright (c) Microsoft Corporation. All rights reserved.
//--------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------
// Constant Buffer Variables
//--------------------------------------------------------------------------------------
Texture2D txDiffuse : register( t0 );
SamplerState samLinear : register( s0 );

cbuffer cbNeverChanges : register( b0 )
{
    matrix View;
};

cbuffer cbChangeOnResize : register( b1 )
{
    matrix Projection;
};

cbuffer cbChangesEveryFrame : register( b2 )
{
    matrix World;
    float4 vMeshColor;
    float4 vLightPos;      // position of the point light source (Lab 5)
    float4 vLightColor;    // light color
    float4 vEyePos;        // camera position
    float4 vAmbient;       // ambient term
    float4 vSpecular;      // rgb - highlight color, w - shininess
    float4 vAttenuation;   // attenuation: 1 / ( x + y * d + z * d * d )
};


//--------------------------------------------------------------------------------------
struct VS_INPUT
{
    float4 Pos : POSITION;
    float3 Norm : NORMAL;
    float2 Tex : TEXCOORD0;
};

struct PS_INPUT
{
    float4 Pos : SV_POSITION;
    float3 Norm : TEXCOORD1;
    float3 PosW : TEXCOORD2;   // world position of the point (Lab 5)
    float2 Tex : TEXCOORD0;
};


//--------------------------------------------------------------------------------------
// Vertex Shader
//--------------------------------------------------------------------------------------
PS_INPUT VS( VS_INPUT input )
{
    PS_INPUT output = (PS_INPUT)0;
    float4 worldPos = mul( input.Pos, World );

    output.Pos = mul( worldPos, View );
    output.Pos = mul( output.Pos, Projection );
    output.PosW = worldPos.xyz;
    // w = 0: the normal rotates with the cube but is not translated
    output.Norm = mul( float4( input.Norm, 0 ), World ).xyz;
    output.Tex = input.Tex;

    return output;
}


//--------------------------------------------------------------------------------------
// Pixel Shader: texture + point light source (Lab 5)
//--------------------------------------------------------------------------------------
float4 PS( PS_INPUT input ) : SV_Target
{
    float4 texColor = txDiffuse.Sample( samLinear, input.Tex ) * vMeshColor;

    float3 N = normalize( input.Norm );                  // normal
    float3 toLight = vLightPos.xyz - input.PosW;         // vector to the light
    float  dist = length( toLight );                     // distance to the light
    float3 L = toLight / max( dist, 1e-5f );             // direction to the light
    float3 V = normalize( vEyePos.xyz - input.PosW );    // direction to the camera
    float3 R = reflect( -L, N );                         // reflected ray

    // distance attenuation of the point light
    float att = 1.0f / ( vAttenuation.x + vAttenuation.y * dist + vAttenuation.z * dist * dist );

    // diffuse term
    float  NdotL = saturate( dot( N, L ) );
    float4 diffuse = NdotL * att * vLightColor;

    // specular term (Phong model), only on the lit side
    float  RdotV = max( 0.0f, dot( R, V ) );
    float4 specular = pow( RdotV, vSpecular.w ) * att * float4( vSpecular.rgb, 0.0f ) * vLightColor;
    specular *= ( NdotL > 0.0f ) ? 1.0f : 0.0f;

    // the texture is modulated by the ambient and diffuse terms, the highlight is added on top
    float4 finalColor = saturate( texColor * ( vAmbient + diffuse ) + specular );
    finalColor.a = 1;
    return finalColor;
}


//--------------------------------------------------------------------------------------
// PSSolid - marker showing the light source position (Lab 5)
//--------------------------------------------------------------------------------------
float4 PSSolid( PS_INPUT input ) : SV_Target
{
    return vMeshColor;
}
