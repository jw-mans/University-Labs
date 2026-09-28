//--------------------------------------------------------------------------------------
// File: Lab4.fx
//
// Lab 4. Author: Andreev Daniil
// Based on: Tutorial06.fx
//
// The cube is lit by a single static light using the Phong model:
// ambient + diffuse (N dot L) + specular (via reflect).
//
// Copyright (c) Microsoft Corporation. All rights reserved.
//--------------------------------------------------------------------------------------


//--------------------------------------------------------------------------------------
// Constant Buffer Variables
//--------------------------------------------------------------------------------------
cbuffer ConstantBuffer : register( b0 )
{
	matrix World;
	matrix View;
	matrix Projection;
	float4 vLightDir;      // direction TO the light source
	float4 vLightColor;    // light color
	float4 vOutputColor;   // color of the light marker
	float4 vEyePos;        // xyz - camera position in world space (Lab 4)
	float4 vAmbient;       // ambient term (Lab 4)
	float4 vDiffuse;       // material color (Lab 4)
	float4 vSpecular;      // rgb - highlight color, w - shininess (Lab 4)
}


//--------------------------------------------------------------------------------------
struct VS_INPUT
{
    float4 Pos : POSITION;
    float3 Norm : NORMAL;
};

struct PS_INPUT
{
    float4 Pos : SV_POSITION;
    float3 Norm : TEXCOORD0;
    float3 PosW : TEXCOORD1;   // world position of the point (Lab 4, needed for the specular term)
};


//--------------------------------------------------------------------------------------
// Vertex Shader
//--------------------------------------------------------------------------------------
PS_INPUT VS( VS_INPUT input )
{
    PS_INPUT output = (PS_INPUT)0;
    output.Pos = mul( input.Pos, World );
    output.PosW = output.Pos.xyz;
    output.Pos = mul( output.Pos, View );
    output.Pos = mul( output.Pos, Projection );
    // w = 0: the normal rotates with the object but is not translated
    output.Norm = mul( float4( input.Norm, 0 ), World ).xyz;

    return output;
}


//--------------------------------------------------------------------------------------
// Pixel Shader: ambient + diffuse + specular (Lab 4)
//--------------------------------------------------------------------------------------
float4 PS( PS_INPUT input ) : SV_Target
{
    float3 N = normalize( input.Norm );                     // normal
    float3 L = normalize( (float3)vLightDir );              // direction to the light
    float3 V = normalize( vEyePos.xyz - input.PosW );       // direction to the camera
    float3 R = reflect( -L, N );                            // reflected ray

    // ambient term
    float4 ambient = vAmbient * vDiffuse;

    // diffuse term (Lambert law)
    float NdotL = saturate( dot( N, L ) );
    float4 diffuse = NdotL * vDiffuse * vLightColor;

    // specular term (Phong model), only on the lit side
    float RdotV = max( 0.0f, dot( R, V ) );
    float4 specular = pow( RdotV, vSpecular.w ) * float4( vSpecular.rgb, 0.0f ) * vLightColor;
    specular *= ( NdotL > 0.0f ) ? 1.0f : 0.0f;

    float4 finalColor = saturate( ambient + diffuse + specular );
    finalColor.a = 1;
    return finalColor;
}


//--------------------------------------------------------------------------------------
// PSSolid - render a solid color (light source marker)
//--------------------------------------------------------------------------------------
float4 PSSolid( PS_INPUT input ) : SV_Target
{
    return vOutputColor;
}
