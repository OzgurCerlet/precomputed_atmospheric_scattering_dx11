/**
 * Copyright (c) 2017 Eric Bruneton
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holders nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */
typedef float scalar;
typedef float3 scalar3;

#include "..\atmosphere_definitions.h"
#include "atmosphere_model.hlsli"
#include "atmosphere_functions.hlsli"

cbuffer atmosphere_cb                   : register(b0)
{
    float4x4 view_from_clip;
    float4x4 world_from_view;
    float3 camera;
    float sun_disk_size_x;
    float3 earth_center;
    float sun_disk_size_y;
    float3 sun_direction;
    float exposure;
    float3 white_point;
};

Texture2D transmittance_texture         : register(t0);
Texture3D scattering_texture            : register(t1);
Texture3D single_mie_scattering_texture : register(t2);
Texture2D irradiance_texture            : register(t3);

SamplerState atmosphere_sampler         : register(s0);

struct PsInput
{
    float4 pos : SV_POSITION;
    float2 uv : TEXCOORD0;
    float3 view_ray_ws : VIEW_RAY_WS;
};

struct PsOutput
{
    float4 color : SV_TARGET;
};

static const float3 kSphereCenter = float3(0.0, 0.0, 1000.0) / kLengthUnitInMeters;
static const float kSphereRadius = 1000.0 / kLengthUnitInMeters;
static const float3 kSphereAlbedo = (float3)0.8;
static const float3 kGroundAlbedo = float3(0.0, 0.0, 0.04);

#ifdef USE_LUMINANCE
#define GetSolarRadiance GetSolarLuminance
#define GetSkyRadiance GetSkyLuminance
#define GetSkyRadianceToPoint GetSkyLuminanceToPoint
#define GetSunAndSkyIrradiance GetSunAndSkyIlluminance
#endif

float GetSunVisibility(float3 pos, float3 sun_direction)
{
    float3 p = pos - kSphereCenter;
    float p_dot_v = dot(p, sun_direction);
    float p_dot_p = dot(p, p);
    float ray_sphere_center_squared_distance = p_dot_p - p_dot_v * p_dot_v;
    float distance_to_intersection = -p_dot_v - sqrt(kSphereRadius * kSphereRadius - ray_sphere_center_squared_distance);
    if (distance_to_intersection > 0.0)
    {
        // Compute the distance between the view ray and the sphere, and the
        // corresponding (tangent of the) subtended angle. Finally, use this to
        // compute an approximate sun visibility.
        float ray_sphere_distance = kSphereRadius - sqrt(ray_sphere_center_squared_distance);
        float ray_sphere_angular_distance = -ray_sphere_distance / p_dot_v;
        return smoothstep(1.0, 0.0, ray_sphere_angular_distance / sun_disk_size_x);
    }
    return 1.0;
}

float GetSkyVisibility(float3 pos)
{
    float3 p = pos - kSphereCenter;
    float p_dot_p = dot(p, p);
    return 1.0 + p.z / sqrt(p_dot_p) * kSphereRadius * kSphereRadius / p_dot_p;
}

void GetSphereShadowInOut(float3 view_direction, float3 sun_direction, out float d_in, out float d_out)
{
    float3 pos = camera - kSphereCenter;
    float pos_dot_sun = dot(pos, sun_direction);
    float view_dot_sun = dot(view_direction, sun_direction);
    float k = sun_disk_size_x;
    float l = 1.0 + k * k;
    float a = 1.0 - l * view_dot_sun * view_dot_sun;
    float b = dot(pos, view_direction) - l * pos_dot_sun * view_dot_sun - k * kSphereRadius * view_dot_sun;
    float c = dot(pos, pos) - l * pos_dot_sun * pos_dot_sun - 2.0 * k * kSphereRadius * pos_dot_sun - kSphereRadius * kSphereRadius;
    float discriminant = b * b - a * c;
    if (discriminant > 0.0)
    {
        d_in = max(0.0, (-b - sqrt(discriminant)) / a);
        d_out = (-b + sqrt(discriminant)) / a;
        // The values of d for which delta is equal to 0 and kSphereRadius / k.
        float d_base = -pos_dot_sun / view_dot_sun;
        float d_apex = -(pos_dot_sun + kSphereRadius / k) / view_dot_sun;
        if (view_dot_sun > 0.0)
        {
            d_in = max(d_in, d_apex);
            d_out = a > 0.0 ? min(d_out, d_base) : d_base;
        }
        else
        {
            d_in = a > 0.0 ? max(d_in, d_base) : d_base;
            d_out = min(d_out, d_apex);
        }
    }
    else
    {
        d_in = 0.0;
        d_out = 0.0;
    }
}

PsOutput ps_main(PsInput input)
{
    #if 1
    // Normalized view direction vector.
    float3 view_direction = normalize(input.view_ray_ws);
    // Tangent of the angle subtended by this fragment.
    float fragment_angular_size = length(ddx(input.view_ray_ws) + ddy(input.view_ray_ws)) / length(input.view_ray_ws);

    float shadow_in;
    float shadow_out;
    GetSphereShadowInOut(view_direction, sun_direction, shadow_in, shadow_out);
    
    // Hack to fade out light shafts when the Sun is very close to the horizon.
    float lightshaft_fadein_hack = smoothstep(0.02, 0.04, dot(normalize(camera - earth_center), sun_direction));

    // Compute the distance between the view ray line and the sphere center,
    // and the distance between the camera and the intersection of the view
    // ray with the sphere (or NaN if there is no intersection).
    float3 p = camera - kSphereCenter;
    float p_dot_v = dot(p, view_direction);
    float p_dot_p = dot(p, p);
    float ray_sphere_center_squared_distance = p_dot_p - p_dot_v * p_dot_v;
    float distance_to_intersection = -p_dot_v - sqrt(kSphereRadius * kSphereRadius - ray_sphere_center_squared_distance);
    
    // Compute the radiance reflected by the sphere, if the ray intersects it.
    float sphere_alpha = 0.0;
    float3 sphere_radiance = (float3) 0.0;
    if (distance_to_intersection > 0.0)
    {
        // Compute the distance between the view ray and the sphere, and the
        // corresponding (tangent of the) subtended angle. Finally, use this to
        // compute the approximate analytic antialiasing factor sphere_alpha.
        float ray_sphere_distance = kSphereRadius - sqrt(ray_sphere_center_squared_distance);
        float ray_sphere_angular_distance = -ray_sphere_distance / p_dot_v;
        sphere_alpha = min(ray_sphere_angular_distance / fragment_angular_size, 1.0);

        float3 pos = camera + view_direction * distance_to_intersection;
        float3 normal = normalize(pos - kSphereCenter);
        
        // Compute the radiance reflected by the sphere.
        float3 sky_irradiance = float3(0, 0, 0);
        float3 sun_irradiance = GetSunAndSkyIrradiance(atmosphere, transmittance_texture, irradiance_texture, atmosphere_sampler, pos - earth_center, normal, sun_direction, sky_irradiance);
        sphere_radiance = kSphereAlbedo * (1.0 / PI) * (sun_irradiance + sky_irradiance);

        float shadow_length = max(0.0, min(shadow_out, distance_to_intersection) - shadow_in) * lightshaft_fadein_hack;
        float3 transmittance = float3(0, 0, 0);
        float3 in_scatter = GetSkyRadianceToPoint(atmosphere, transmittance_texture, scattering_texture, single_mie_scattering_texture, atmosphere_sampler, camera - earth_center, pos - earth_center, shadow_length, sun_direction, transmittance);
        sphere_radiance = sphere_radiance * transmittance + in_scatter;
    }
    
    // Compute the distance between the view ray line and the Earth center,
    // and the distance between the camera and the intersection of the view
    // ray with the ground (or NaN if there is no intersection).
    p = camera - earth_center;
    p_dot_v = dot(p, view_direction);
    p_dot_p = dot(p, p);
    float ray_earth_center_squared_distance = p_dot_p - p_dot_v * p_dot_v;
    distance_to_intersection = -p_dot_v - sqrt(earth_center.z * earth_center.z - ray_earth_center_squared_distance);
    
    // Compute the radiance reflected by the ground, if the ray intersects it.
    float ground_alpha = 0.0;
    float3 ground_radiance = (float3) 0.0;
    if (distance_to_intersection > 0.0)
    {
        float3 pos = camera + view_direction * distance_to_intersection;
        float3 normal = normalize(pos - earth_center);
        
        // Compute the radiance reflected by the ground.
        float3 sky_irradiance = float3(0,0,0);
        float3 sun_irradiance = GetSunAndSkyIrradiance(atmosphere, transmittance_texture, irradiance_texture, atmosphere_sampler, pos - earth_center, normal, sun_direction, sky_irradiance);
        ground_radiance = kGroundAlbedo * (1.0 / PI) * (sun_irradiance * GetSunVisibility(pos, sun_direction) + sky_irradiance * GetSkyVisibility(pos));

        float shadow_length = max(0.0, min(shadow_out, distance_to_intersection) - shadow_in) * lightshaft_fadein_hack;
        float3 transmittance = float3(0, 0, 0);
        float3 in_scatter = GetSkyRadianceToPoint(atmosphere, transmittance_texture, scattering_texture, single_mie_scattering_texture, atmosphere_sampler, camera - earth_center, pos - earth_center, shadow_length, sun_direction, transmittance);
        ground_radiance = ground_radiance * transmittance + in_scatter;
        ground_alpha = 1.0;
    }
    
    // Compute the radiance of the sky.
    float shadow_length = max(0.0, shadow_out - shadow_in) * lightshaft_fadein_hack;
    float3 transmittance;
    float3 radiance = GetSkyRadiance(atmosphere, transmittance_texture, scattering_texture, single_mie_scattering_texture, atmosphere_sampler, camera - earth_center, view_direction, shadow_length, sun_direction, transmittance);
    
    // If the view ray intersects the Sun, add the Sun radiance.
    if (dot(view_direction, sun_direction) > sun_disk_size_y) {
        radiance = radiance + transmittance * GetSolarRadiance(atmosphere);
    }
    
    radiance = lerp(radiance, ground_radiance, ground_alpha);
    radiance = lerp(radiance, sphere_radiance, sphere_alpha);

    float3 color = pow((float3) 1.0 - exp(-radiance / white_point * exposure), (float3)(1.0 / 2.2));

#else

    //float3 view_direction = normalize(input.view_ray_ws);
    float3 view_direction = float3(-0.34738, 0.00081, -0.93772);
    float shadow_length = 0;
    float3 transmittance;
    //float3 radiance = GetSkyRadiance(atmosphere, transmittance_texture, scattering_texture, single_mie_scattering_texture, atmosphere_sampler, float3(9.00, 0.00, 6359.99805), float3(-0.99829, -0.03879, -0.04376), 0.0, float3(0.00059, 0.0, -1.0), transmittance);
    float3 radiance = GetSkyRadiance(atmosphere, transmittance_texture, scattering_texture, single_mie_scattering_texture, atmosphere_sampler, camera - earth_center, view_direction, shadow_length, sun_direction, transmittance);
    
    float3 color = radiance;
    //float3 color = pow((float3) 1.0 - exp(-radiance / white_point * exposure), (float3) (1.0 / 2.2));

#endif

    PsOutput output = (PsOutput) 0;
    output.color.xyz = color;
    return output;
}