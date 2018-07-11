typedef float scalar;
typedef float3 scalar3;

#include "..\atmosphere_definitions.h"
#include "atmosphere_model.hlsli"
#include "atmosphere_functions.hlsli"

#pragma pack_matrix(row_major)
cbuffer precompute_cb : register(b0)
{
    float3x3 luminance_from_radiance : packoffset(c0);
    int scattering_order : packoffset(c4.x);
};

struct PsInput {
#if(PRECOMPUTE_SHADER_TYPE==TRANSMITTANCE)
    float4  pos_ss      : SV_POSITION;
    float2  uv          : TEXCOORD0;
#elif(PRECOMPUTE_SHADER_TYPE==DELTA_IRRADIANCE)
    float4  pos_ss      : SV_POSITION;
    float2  uv          : TEXCOORD0;
#elif(PRECOMPUTE_SHADER_TYPE==SINGLE_SCATTERING)
    float4  pos_ss      : SV_POSITION;
    float2  uv          : TEXCOORD0;
    uint    slice_index : SV_RENDERTARGETARRAYINDEX;
#elif(PRECOMPUTE_SHADER_TYPE==SCATTERING_DENSITY)
    float4  pos_ss      : SV_POSITION;
    float2  uv          : TEXCOORD0;
    uint    slice_index : SV_RENDERTARGETARRAYINDEX;
#elif(PRECOMPUTE_SHADER_TYPE==INDIRECT_IRRADIANCE)
    float4  pos_ss      : SV_POSITION;
    float2  uv          : TEXCOORD0;
    uint    slice_index : SV_RENDERTARGETARRAYINDEX;
#elif(PRECOMPUTE_SHADER_TYPE==MULTIPLE_SCATTERING)
    float4  pos_ss      : SV_POSITION;
    float2  uv          : TEXCOORD0;
    uint    slice_index : SV_RENDERTARGETARRAYINDEX;
#endif
};

struct PsOutput
{
#if(PRECOMPUTE_SHADER_TYPE==TRANSMITTANCE)
    float4 transmittance                : SV_TARGET0;
#elif(PRECOMPUTE_SHADER_TYPE==DELTA_IRRADIANCE)
    float3 delta_irradiance             : SV_TARGET0;
    float3 irradiance                   : SV_TARGET1;
#elif(PRECOMPUTE_SHADER_TYPE==SINGLE_SCATTERING)
    float3 delta_rayleigh               : SV_TARGET0;
    float3 delta_mie                    : SV_TARGET1;
    float4 scattering                   : SV_TARGET2;
    float3 single_mie                   : SV_TARGET3;
#elif(PRECOMPUTE_SHADER_TYPE==SCATTERING_DENSITY)
    float3 scattering_density           : SV_TARGET0;
#elif(PRECOMPUTE_SHADER_TYPE==INDIRECT_IRRADIANCE)
    float3 delta_irradiance             : SV_TARGET0;
    float3 irradiance                   : SV_TARGET1;
#elif(PRECOMPUTE_SHADER_TYPE==MULTIPLE_SCATTERING)
    float3 delta_multiple_scattering    : SV_TARGET0;
    float4 scattering                   : SV_TARGET1;
#endif
};

#if(PRECOMPUTE_SHADER_TYPE==TRANSMITTANCE)
Texture2D transmittance_texture                 : register(t0);
#elif(PRECOMPUTE_SHADER_TYPE==DELTA_IRRADIANCE)
Texture2D transmittance_texture                 : register(t0);
#elif(PRECOMPUTE_SHADER_TYPE==SINGLE_SCATTERING)
Texture2D transmittance_texture                 : register(t0);
#elif(PRECOMPUTE_SHADER_TYPE==SCATTERING_DENSITY)
Texture2D transmittance_texture                 : register(t0);
Texture3D single_rayleigh_scattering_texture    : register(t1);
Texture3D single_mie_scattering_texture         : register(t2);
Texture3D multiple_scattering_texture           : register(t3);
Texture2D delta_irradiance_texture              : register(t4);
#elif(PRECOMPUTE_SHADER_TYPE==INDIRECT_IRRADIANCE)
Texture3D single_rayleigh_scattering_texture    : register(t0);
Texture3D single_mie_scattering_texture         : register(t1);
Texture3D multiple_scattering_texture           : register(t2);
#elif(PRECOMPUTE_SHADER_TYPE==MULTIPLE_SCATTERING)
Texture2D transmittance_texture                 : register(t0);
Texture3D scattering_density_texture            : register(t1);
#endif

SamplerState atmosphere_sampler : register(s0);

PsOutput ps_main(PsInput input) {
    PsOutput output = (PsOutput) 0;
    float2 coord_ss = input.pos_ss.xy;
#if(PRECOMPUTE_SHADER_TYPE==TRANSMITTANCE)
    output.transmittance = float4(ComputeTransmittanceToTopAtmosphereBoundaryTexture(atmosphere, coord_ss), 1.0);
#elif(PRECOMPUTE_SHADER_TYPE==DELTA_IRRADIANCE)
    output.delta_irradiance = ComputeDirectIrradianceTexture(atmosphere,transmittance_texture, atmosphere_sampler, coord_ss);
    output.irradiance = (float3) 0.0;
#elif(PRECOMPUTE_SHADER_TYPE==SINGLE_SCATTERING)
    ComputeSingleScatteringTexture(atmosphere, transmittance_texture, atmosphere_sampler, float3(coord_ss, input.slice_index + 0.5), output.delta_rayleigh, output.delta_mie);
    output.scattering = float4(mul(luminance_from_radiance, output.delta_rayleigh), mul(luminance_from_radiance, output.delta_mie).r);
    output.single_mie = mul(luminance_from_radiance, output.delta_mie);
#elif(PRECOMPUTE_SHADER_TYPE==SCATTERING_DENSITY)
    output.scattering_density = ComputeScatteringDensityTexture(atmosphere, transmittance_texture, single_rayleigh_scattering_texture, single_mie_scattering_texture, 
    multiple_scattering_texture, delta_irradiance_texture, atmosphere_sampler, float3(coord_ss, input.slice_index + 0.5), scattering_order);
#elif(PRECOMPUTE_SHADER_TYPE==INDIRECT_IRRADIANCE)
    output.delta_irradiance = ComputeIndirectIrradianceTexture(atmosphere, single_rayleigh_scattering_texture, single_mie_scattering_texture, multiple_scattering_texture, atmosphere_sampler, float3(coord_ss, input.slice_index + 0.5), scattering_order);
    output.irradiance = mul(luminance_from_radiance, output.delta_irradiance);
#elif(PRECOMPUTE_SHADER_TYPE==MULTIPLE_SCATTERING)
    float nu = 0.0;
    output.delta_multiple_scattering = ComputeMultipleScatteringTexture(atmosphere, transmittance_texture, scattering_density_texture, atmosphere_sampler, float3(coord_ss, input.slice_index + 0.5), nu);
    output.scattering = float4(mul(luminance_from_radiance, output.delta_multiple_scattering / RayleighPhaseFunction(nu)), 0.0);
#endif
    return output;
}