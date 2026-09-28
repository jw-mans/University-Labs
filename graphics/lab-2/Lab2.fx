//--------------------------------------------------------------------------------------
// File: Lab2.fx
//
// Lab 2. Author: Andreev Daniil
// Based on: Tutorial02.fx
//--------------------------------------------------------------------------------------

// Translation and scale of the quad (Lab 2 assignment) - applied in the vertex shader
static const float3 QUAD_SCALE  = float3( 1.4f, 1.4f, 1.0f );
static const float3 QUAD_OFFSET = float3( 0.25f, 0.15f, 0.0f );

// Quad color (Tutorial02 used yellow - float4(1, 1, 0, 1))
static const float4 QUAD_COLOR = float4( 0.95f, 0.45f, 0.1f, 1.0f );   // orange

//--------------------------------------------------------------------------------------
// Vertex Shader
//--------------------------------------------------------------------------------------
float4 VS( float4 Pos : POSITION ) : SV_POSITION
{
    float3 p = Pos.xyz * QUAD_SCALE + QUAD_OFFSET;
    return float4( p, 1.0f );
}


//--------------------------------------------------------------------------------------
// Pixel Shader
//--------------------------------------------------------------------------------------
float4 PS( float4 Pos : SV_POSITION ) : SV_Target
{
    return QUAD_COLOR;
}
