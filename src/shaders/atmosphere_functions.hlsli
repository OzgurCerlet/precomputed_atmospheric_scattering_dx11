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
 *
 * Precomputed Atmospheric Scattering
 * Copyright (c) 2008 INRIA
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

float mod(float x, float y)
{
    return x - y * floor(x / y);
}

float ClampCosine(float mu)
{
    return clamp(mu, float(-1.0), float(1.0));
}

float ClampDistance(float d)
{
    return max(d, 0.0);
}

float ClampRadius(const in AtmosphereParameters atmosphere, float r)
{
    return clamp(r, atmosphere.bottom_radius, atmosphere.top_radius);
}

float SafeSqrt(float a)
{
    return sqrt(max(a, 0.0));
}

float DistanceToTopAtmosphereBoundary(const in AtmosphereParameters atmosphere, float r, float mu)
{
    float discriminant = r * r * (mu * mu - 1.0) + atmosphere.top_radius * atmosphere.top_radius;
    return ClampDistance(-r * mu + SafeSqrt(discriminant));
}

float DistanceToBottomAtmosphereBoundary(const in AtmosphereParameters atmosphere, float r, float mu)
{
    float discriminant = r * r * (mu * mu - 1.0) + atmosphere.bottom_radius * atmosphere.bottom_radius;
    return ClampDistance(-r * mu - SafeSqrt(discriminant));
}

bool RayIntersectsGround(const in AtmosphereParameters atmosphere, float r, float mu)
{
    return (mu < 0.0) && ((r * r * (mu * mu - 1.0) + atmosphere.bottom_radius * atmosphere.bottom_radius) >= 0.0);
}

float GetLayerDensity(const in DensityProfileLayer layer, float altitude)
{
    float density = layer.exp_term * exp(layer.exp_scale * altitude) + layer.linear_term * altitude + layer.constant_term;
    return clamp(density, float(0.0), float(1.0));
}

float GetProfileDensity(const in DensityProfile profile, float altitude)
{
    return (altitude < profile.layers[0].width) ? GetLayerDensity(profile.layers[0], altitude) : GetLayerDensity(profile.layers[1], altitude);
}

float ComputeOpticalLengthToTopAtmosphereBoundary(const in AtmosphereParameters atmosphere, const in DensityProfile profile, float r, float mu)
{
    const int SAMPLE_COUNT = 500;
    float dx = DistanceToTopAtmosphereBoundary(atmosphere, r, mu) / float(SAMPLE_COUNT);
    float result = 0.0;
    for (int i = 0; i <= SAMPLE_COUNT; ++i)
    {
        float d_i = float(i) * dx;
        float r_i = sqrt(d_i * d_i + 2.0 * r * mu * d_i + r * r);
        float y_i = GetProfileDensity(profile, r_i - atmosphere.bottom_radius);
        float weight_i = (i == 0 || i == SAMPLE_COUNT) ? 0.5 : 1.0;
        result += y_i * weight_i * dx;
    }
    return result;
}

float3 ComputeTransmittanceToTopAtmosphereBoundary(const in AtmosphereParameters atmosphere, float r, float mu)
{
    return exp(-(
			atmosphere.rayleigh_scattering *
			ComputeOpticalLengthToTopAtmosphereBoundary(atmosphere, atmosphere.rayleigh_density, r, mu) + atmosphere.mie_extinction *
			ComputeOpticalLengthToTopAtmosphereBoundary(atmosphere, atmosphere.mie_density, r, mu) + atmosphere.absorption_extinction *
			ComputeOpticalLengthToTopAtmosphereBoundary(atmosphere, atmosphere.absorption_density, r, mu))
	);
}

float GetTextureCoordFromUnitRange(float x, int texture_size)
{
    return 0.5 / float(texture_size) + x * (1.0 - 1.0 / float(texture_size));
}

float GetUnitRangeFromTextureCoord(float u, int texture_size)
{
    return (u - 0.5 / float(texture_size)) / (1.0 - 1.0 / float(texture_size));
}

float2 GetTransmittanceTextureUvFromRMu(const in AtmosphereParameters atmosphere, float r, float mu)
{
    float H = sqrt(atmosphere.top_radius * atmosphere.top_radius - atmosphere.bottom_radius * atmosphere.bottom_radius);
    float rho = SafeSqrt(r * r - atmosphere.bottom_radius * atmosphere.bottom_radius);
    float d = DistanceToTopAtmosphereBoundary(atmosphere, r, mu);
    float d_min = atmosphere.top_radius - r;
    float d_max = rho + H;
    float x_mu = (d - d_min) / (d_max - d_min);
    float x_r = rho / H;
    return float2(GetTextureCoordFromUnitRange(x_mu, TRANSMITTANCE_TEXTURE_WIDTH), GetTextureCoordFromUnitRange(x_r, TRANSMITTANCE_TEXTURE_HEIGHT));
}

void GetRMuFromTransmittanceTextureUv(const in AtmosphereParameters atmosphere, const in float2 uv, out float r, out float mu)
{
    float x_mu = GetUnitRangeFromTextureCoord(uv.x, TRANSMITTANCE_TEXTURE_WIDTH);
    float x_r = GetUnitRangeFromTextureCoord(uv.y, TRANSMITTANCE_TEXTURE_HEIGHT);
    float H = sqrt(atmosphere.top_radius * atmosphere.top_radius - atmosphere.bottom_radius * atmosphere.bottom_radius);
    float rho = H * x_r;
    r = sqrt(rho * rho + atmosphere.bottom_radius * atmosphere.bottom_radius);
    float d_min = atmosphere.top_radius - r;
    float d_max = rho + H;
    float d = d_min + x_mu * (d_max - d_min);
    mu = d == 0.0 ? float(1.0) : (H * H - rho * rho - d * d) / (2.0 * r * d);
    mu = ClampCosine(mu);
}

float3 ComputeTransmittanceToTopAtmosphereBoundaryTexture(const in AtmosphereParameters atmosphere, const in float2 frag_coord)
{
    const float2 TRANSMITTANCE_TEXTURE_SIZE = float2(TRANSMITTANCE_TEXTURE_WIDTH, TRANSMITTANCE_TEXTURE_HEIGHT);
    float r;
    float mu;
    GetRMuFromTransmittanceTextureUv(atmosphere, frag_coord / TRANSMITTANCE_TEXTURE_SIZE, r, mu);
    return ComputeTransmittanceToTopAtmosphereBoundary(atmosphere, r, mu);
}

float3 GetTransmittanceToTopAtmosphereBoundary(const in AtmosphereParameters atmosphere, const in Texture2D transmittance_texture, const in SamplerState atmo_sampler, float r, float mu)
{
    float2 uv = GetTransmittanceTextureUvFromRMu(atmosphere, r, mu);
    return transmittance_texture.SampleLevel(atmo_sampler, uv, 0).xyz;
}

float3 GetTransmittance(const in AtmosphereParameters atmosphere, const in Texture2D transmittance_texture, const in SamplerState atmo_sampler, float r, float mu, float d, bool ray_r_mu_intersects_ground)
{
    float r_d = ClampRadius(atmosphere, sqrt(d * d + 2.0 * r * mu * d + r * r));
    float mu_d = ClampCosine((r * mu + d) / r_d);
    if (ray_r_mu_intersects_ground)
    {
        return min(GetTransmittanceToTopAtmosphereBoundary(atmosphere, transmittance_texture, atmo_sampler, r_d, -mu_d) / GetTransmittanceToTopAtmosphereBoundary(atmosphere, transmittance_texture, atmo_sampler, r, -mu), (float3) (1.0));
    }
    else
    {
        return min(GetTransmittanceToTopAtmosphereBoundary(atmosphere, transmittance_texture, atmo_sampler, r, mu) / GetTransmittanceToTopAtmosphereBoundary(atmosphere, transmittance_texture, atmo_sampler, r_d, mu_d), (float3) (1.0));
    }
}

float3 GetTransmittanceToSun(const in AtmosphereParameters atmosphere, const in Texture2D transmittance_texture, const in SamplerState atmo_sampler, float r, float mu_s)
{
    float sin_theta_h = atmosphere.bottom_radius / r;
    float cos_theta_h = -sqrt(max(1.0 - sin_theta_h * sin_theta_h, 0.0));
    return GetTransmittanceToTopAtmosphereBoundary(atmosphere, transmittance_texture, atmo_sampler, r, mu_s) * smoothstep(-sin_theta_h * atmosphere.sun_angular_radius, sin_theta_h * atmosphere.sun_angular_radius, mu_s - cos_theta_h);
}

void ComputeSingleScatteringIntegrand(
	const in AtmosphereParameters atmosphere,
	const in Texture2D transmittance_texture, SamplerState atmo_sampler,
	float r, float mu, float mu_s, float nu, float d,
	bool ray_r_mu_intersects_ground,
	out float3 rayleigh, out float3 mie)
{
    float r_d = ClampRadius(atmosphere, sqrt(d * d + 2.0 * r * mu * d + r * r));
    float mu_s_d = ClampCosine((r * mu_s + d * nu) / r_d);
    float3 transmittance = GetTransmittance(atmosphere, transmittance_texture, atmo_sampler, r, mu, d, ray_r_mu_intersects_ground) * GetTransmittanceToSun(atmosphere, transmittance_texture, atmo_sampler, r_d, mu_s_d);
    rayleigh = transmittance * GetProfileDensity(atmosphere.rayleigh_density, r_d - atmosphere.bottom_radius);
    mie = transmittance * GetProfileDensity(atmosphere.mie_density, r_d - atmosphere.bottom_radius);
}

float DistanceToNearestAtmosphereBoundary(const in AtmosphereParameters atmosphere, float r, float mu, bool ray_r_mu_intersects_ground)
{
    if (ray_r_mu_intersects_ground)
    {
        return DistanceToBottomAtmosphereBoundary(atmosphere, r, mu);
    }
    else
    {
        return DistanceToTopAtmosphereBoundary(atmosphere, r, mu);
    }
}

void ComputeSingleScattering(
	const in AtmosphereParameters atmosphere,
	const in Texture2D transmittance_texture, const in SamplerState atmo_sampler,
	float r, float mu, float mu_s, float nu,
	bool ray_r_mu_intersects_ground,
	out float3 rayleigh, out float3 mie)
{
    const int SAMPLE_COUNT = 50;
    float dx = DistanceToNearestAtmosphereBoundary(atmosphere, r, mu, ray_r_mu_intersects_ground) / float(SAMPLE_COUNT);
    float3 rayleigh_sum = (float3) (0.0);
    float3 mie_sum = (float3) (0.0);
    for (int i = 0; i <= SAMPLE_COUNT; ++i)
    {
        float d_i = float(i) * dx;
        float3 rayleigh_i;
        float3 mie_i;
        ComputeSingleScatteringIntegrand(atmosphere, transmittance_texture, atmo_sampler, r, mu, mu_s, nu, d_i, ray_r_mu_intersects_ground, rayleigh_i, mie_i);
        float weight_i = (i == 0 || i == SAMPLE_COUNT) ? 0.5 : 1.0;
        rayleigh_sum += rayleigh_i * weight_i;
        mie_sum += mie_i * weight_i;
    }
    rayleigh = rayleigh_sum * dx * atmosphere.solar_irradiance * atmosphere.rayleigh_scattering;
    mie = mie_sum * dx * atmosphere.solar_irradiance * atmosphere.mie_scattering;
}

float RayleighPhaseFunction(float nu)
{
    float k = 3.0 / (16.0 * PI);
    return k * (1.0 + nu * nu);
}

float MiePhaseFunction(float g, float nu)
{
    float k = 3.0 / (8.0 * PI) * (1.0 - g * g) / (2.0 + g * g);
    return k * (1.0 + nu * nu) / pow(1.0 + g * g - 2.0 * g * nu, 1.5);
}

float4 GetScatteringTextureUvwzFromRMuMuSNu(
	const in AtmosphereParameters atmosphere,
	float r, float mu, float mu_s, float nu,
	bool ray_r_mu_intersects_ground)
{
    float H = sqrt(atmosphere.top_radius * atmosphere.top_radius - atmosphere.bottom_radius * atmosphere.bottom_radius);
    float rho = SafeSqrt(r * r - atmosphere.bottom_radius * atmosphere.bottom_radius);
    float u_r = GetTextureCoordFromUnitRange(rho / H, SCATTERING_TEXTURE_R_SIZE);
    float r_mu = r * mu;
    float discriminant = r_mu * r_mu - r * r + atmosphere.bottom_radius * atmosphere.bottom_radius;
    float u_mu;
    if (ray_r_mu_intersects_ground)
    {
        float d = -r_mu - SafeSqrt(discriminant);
        float d_min = r - atmosphere.bottom_radius;
        float d_max = rho;
        u_mu = 0.5 - 0.5 * GetTextureCoordFromUnitRange(d_max == d_min ? 0.0 : (d - d_min) / (d_max - d_min), SCATTERING_TEXTURE_MU_SIZE / 2);
    }
    else
    {
        float d = -r_mu + SafeSqrt(discriminant + H * H);
        float d_min = atmosphere.top_radius - r;
        float d_max = rho + H;
        u_mu = 0.5 + 0.5 * GetTextureCoordFromUnitRange((d - d_min) / (d_max - d_min), SCATTERING_TEXTURE_MU_SIZE / 2);
    }

    float d = DistanceToTopAtmosphereBoundary(atmosphere, atmosphere.bottom_radius, mu_s);
    float d_min = atmosphere.top_radius - atmosphere.bottom_radius;
    float d_max = H;
    float a = (d - d_min) / (d_max - d_min);
    float A = -2.0 * atmosphere.mu_s_min * atmosphere.bottom_radius / (d_max - d_min);
    float u_mu_s = GetTextureCoordFromUnitRange(max(1.0 - a / A, 0.0) / (1.0 + a), SCATTERING_TEXTURE_MU_S_SIZE);
    float u_nu = (nu + 1.0) / 2.0;

    return float4(u_nu, u_mu_s, u_mu, u_r);
}

void GetRMuMuSNuFromScatteringTextureUvwz(const in AtmosphereParameters atmosphere,
	const in float4 uvwz, out float r, out float mu, out float mu_s, out float nu, out bool ray_r_mu_intersects_ground)
{
    float H = sqrt(atmosphere.top_radius * atmosphere.top_radius - atmosphere.bottom_radius * atmosphere.bottom_radius);
    float rho = H * GetUnitRangeFromTextureCoord(uvwz.w, SCATTERING_TEXTURE_R_SIZE);
    r = sqrt(rho * rho + atmosphere.bottom_radius * atmosphere.bottom_radius);
    if (uvwz.z < 0.5)
    {
        float d_min = r - atmosphere.bottom_radius;
        float d_max = rho;
        float d = d_min + (d_max - d_min) * GetUnitRangeFromTextureCoord(1.0 - 2.0 * uvwz.z, SCATTERING_TEXTURE_MU_SIZE / 2);
        mu = d == 0.0 ? float(-1.0) : ClampCosine(-(rho * rho + d * d) / (2.0 * r * d));
        ray_r_mu_intersects_ground = true;
    }
    else
    {
        float d_min = atmosphere.top_radius - r;
        float d_max = rho + H;
        float d = d_min + (d_max - d_min) * GetUnitRangeFromTextureCoord(2.0 * uvwz.z - 1.0, SCATTERING_TEXTURE_MU_SIZE / 2);
        mu = d == 0.0 ? float(1.0) : ClampCosine((H * H - rho * rho - d * d) / (2.0 * r * d));
        ray_r_mu_intersects_ground = false;
    }

    float x_mu_s = GetUnitRangeFromTextureCoord(uvwz.y, SCATTERING_TEXTURE_MU_S_SIZE);
    float d_min = atmosphere.top_radius - atmosphere.bottom_radius;
    float d_max = H;
    float A = -2.0 * atmosphere.mu_s_min * atmosphere.bottom_radius / (d_max - d_min);
    float a = (A - x_mu_s * A) / (1.0 + x_mu_s * A);
    float d = d_min + min(a, A) * (d_max - d_min);
    mu_s = d == 0.0 ? float(1.0) : ClampCosine((H * H - d * d) / (2.0 * atmosphere.bottom_radius * d));
    nu = ClampCosine(uvwz.x * 2.0 - 1.0);
}

void GetRMuMuSNuFromScatteringTextureFragCoord(
	const in AtmosphereParameters atmosphere, const in float3 frag_coord,
	out float r, out float mu, out float mu_s, out float nu,
	out bool ray_r_mu_intersects_ground)
{
    const float4 SCATTERING_TEXTURE_SIZE = float4(SCATTERING_TEXTURE_NU_SIZE - 1, SCATTERING_TEXTURE_MU_S_SIZE, SCATTERING_TEXTURE_MU_SIZE, SCATTERING_TEXTURE_R_SIZE);
    float frag_coord_nu = floor(frag_coord.x / float(SCATTERING_TEXTURE_MU_S_SIZE));
    float frag_coord_mu_s = mod(frag_coord.x, float(SCATTERING_TEXTURE_MU_S_SIZE));
    float4 uvwz = float4(frag_coord_nu, frag_coord_mu_s, frag_coord.y, frag_coord.z) / SCATTERING_TEXTURE_SIZE;
    GetRMuMuSNuFromScatteringTextureUvwz(atmosphere, uvwz, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    nu = clamp(nu, mu * mu_s - sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)), mu * mu_s + sqrt((1.0 - mu * mu) * (1.0 - mu_s * mu_s)));
}

void ComputeSingleScatteringTexture(const in AtmosphereParameters atmosphere,
	const in Texture2D transmittance_texture, const in SamplerState atmo_sampler,
	const in float3 frag_coord,
	out float3 rayleigh, out float3 mie)
{
    float r;
    float mu;
    float mu_s;
    float nu;
    bool ray_r_mu_intersects_ground;
    GetRMuMuSNuFromScatteringTextureFragCoord(atmosphere, frag_coord, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    ComputeSingleScattering(atmosphere, transmittance_texture, atmo_sampler, r, mu, mu_s, nu, ray_r_mu_intersects_ground, rayleigh, mie);
}

float3 GetScattering(
	const in AtmosphereParameters atmosphere,
	const in Texture3D scattering_texture, const in SamplerState atmo_sampler,
	float r, float mu, float mu_s, float nu,
	bool ray_r_mu_intersects_ground)
{
    float4 uvwz = GetScatteringTextureUvwzFromRMuMuSNu(atmosphere, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    float tex_coord_x = uvwz.x * float(SCATTERING_TEXTURE_NU_SIZE - 1);
    float tex_x = floor(tex_coord_x);
    float lerp = tex_coord_x - tex_x;
    float3 uvw0 = float3((tex_x + uvwz.y) / float(SCATTERING_TEXTURE_NU_SIZE), uvwz.z, uvwz.w);
    float3 uvw1 = float3((tex_x + 1.0 + uvwz.y) / float(SCATTERING_TEXTURE_NU_SIZE), uvwz.z, uvwz.w);
    return (scattering_texture.Sample(atmo_sampler, uvw0) * (1.0 - lerp) + scattering_texture.Sample(atmo_sampler, uvw1) * lerp).xyz;
}

float3 GetScattering(
	const in AtmosphereParameters atmosphere,
	const in Texture3D single_rayleigh_scattering_texture,
	const in Texture3D single_mie_scattering_texture,
	const in Texture3D multiple_scattering_texture,
	SamplerState atmo_sampler,
	float r, float mu, float mu_s, float nu,
	bool ray_r_mu_intersects_ground,
	int scattering_order)
{
    if (scattering_order == 1)
    {
        float3 rayleigh = GetScattering(atmosphere, single_rayleigh_scattering_texture, atmo_sampler, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
        float3 mie = GetScattering(atmosphere, single_mie_scattering_texture, atmo_sampler, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
        return rayleigh * RayleighPhaseFunction(nu) + mie * MiePhaseFunction(atmosphere.mie_phase_function_g, nu);
    }
    else
    {
        return GetScattering(atmosphere, multiple_scattering_texture, atmo_sampler, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    }
}

float2 GetIrradianceTextureUvFromRMuS(const in AtmosphereParameters atmosphere, float r, float mu_s)
{
    float x_r = (r - atmosphere.bottom_radius) / (atmosphere.top_radius - atmosphere.bottom_radius);
    float x_mu_s = mu_s * 0.5 + 0.5;
    return float2(GetTextureCoordFromUnitRange(x_mu_s, IRRADIANCE_TEXTURE_WIDTH), GetTextureCoordFromUnitRange(x_r, IRRADIANCE_TEXTURE_HEIGHT));
}

float3 GetIrradiance(
	const in AtmosphereParameters atmosphere,
	Texture2D irradiance_texture,
	SamplerState atmo_sampler,
	float r, float mu_s)
{
    float2 uv = GetIrradianceTextureUvFromRMuS(atmosphere, r, mu_s);
    return irradiance_texture.SampleLevel(atmo_sampler, uv, 0).xyz;
}

static const float2 IRRADIANCE_TEXTURE_SIZE = float2(IRRADIANCE_TEXTURE_WIDTH, IRRADIANCE_TEXTURE_HEIGHT);

void GetRMuSFromIrradianceTextureUv(const in AtmosphereParameters atmosphere, float2 uv, out float r, out float mu_s)
{
    float x_mu_s = GetUnitRangeFromTextureCoord(uv.x, IRRADIANCE_TEXTURE_WIDTH);
    float x_r = GetUnitRangeFromTextureCoord(uv.y, IRRADIANCE_TEXTURE_HEIGHT);
    r = atmosphere.bottom_radius + x_r * (atmosphere.top_radius - atmosphere.bottom_radius);
    mu_s = ClampCosine(2.0 * x_mu_s - 1.0);
}

float3 ComputeScatteringDensity(
	const in AtmosphereParameters atmosphere,
	Texture2D transmittance_texture,
	Texture3D single_rayleigh_scattering_texture,
	Texture3D single_mie_scattering_texture,
	Texture3D multiple_scattering_texture,
	Texture2D irradiance_texture,
	SamplerState atmo_sampler,
	float r, float mu, float mu_s, float nu, int scattering_order)
{
    float3 zenith_direction = float3(0.0, 0.0, 1.0);
    float3 omega = float3(sqrt(1.0 - mu * mu), 0.0, mu);
    float sun_dir_x = (omega.x == 0.0) ? 0.0 : (nu - mu * mu_s) / omega.x;
    float sun_dir_y = sqrt(max(1.0 - sun_dir_x * sun_dir_x - mu_s * mu_s, 0.0));
    float3 omega_s = float3(sun_dir_x, sun_dir_y, mu_s);
    const int SAMPLE_COUNT = 16;
    const float dphi = PI / float(SAMPLE_COUNT);
    const float dtheta = PI / float(SAMPLE_COUNT);
    float3 rayleigh_mie = (float3) (0.0);
    for (int l = 0; l < SAMPLE_COUNT; ++l)
    {
        float theta = (float(l) + 0.5) * dtheta;
        float cos_theta = cos(theta);
        float sin_theta = sin(theta);
        bool ray_r_theta_intersects_ground = RayIntersectsGround(atmosphere, r, cos_theta);
        float distance_to_ground = 0.0;
        float3 transmittance_to_ground = (float3) (0.0);
        float3 ground_albedo = (float3) (0.0);
		
        if (ray_r_theta_intersects_ground)
        {
            distance_to_ground = DistanceToBottomAtmosphereBoundary(atmosphere, r, cos_theta);
            transmittance_to_ground = GetTransmittance(atmosphere, transmittance_texture, atmo_sampler, r, cos_theta, distance_to_ground, true /* ray_intersects_ground */);
            ground_albedo = atmosphere.ground_albedo;
        }

        for (int m = 0; m < 2 * SAMPLE_COUNT; ++m)
        {
            float phi = (float(m) + 0.5) * dphi;
            float3 omega_i = float3(cos(phi) * sin_theta, sin(phi) * sin_theta, cos_theta);
            float domega_i = (dtheta) * (dphi) * sin(theta);
            float nu1 = dot(omega_s, omega_i);
            float3 incident_radiance = GetScattering(atmosphere, single_rayleigh_scattering_texture, single_mie_scattering_texture, multiple_scattering_texture, atmo_sampler, r, omega_i.z, mu_s, nu1, ray_r_theta_intersects_ground, scattering_order - 1);
            float3 ground_normal = normalize(zenith_direction * r + omega_i * distance_to_ground);
            float3 ground_irradiance = GetIrradiance(atmosphere, irradiance_texture, atmo_sampler, atmosphere.bottom_radius, dot(ground_normal, omega_s));
            incident_radiance += transmittance_to_ground * ground_albedo * (1.0 / PI) * ground_irradiance;
            float nu2 = dot(omega, omega_i);
            float rayleigh_density = GetProfileDensity(atmosphere.rayleigh_density, r - atmosphere.bottom_radius);
            float mie_density = GetProfileDensity(atmosphere.mie_density, r - atmosphere.bottom_radius);
            rayleigh_mie += incident_radiance * (atmosphere.rayleigh_scattering * rayleigh_density * RayleighPhaseFunction(nu2) + atmosphere.mie_scattering * mie_density * MiePhaseFunction(atmosphere.mie_phase_function_g, nu2)) * domega_i;
        }
    }
    return rayleigh_mie;
}

float3 ComputeMultipleScattering(
	const in AtmosphereParameters atmosphere,
	Texture2D transmittance_texture,
	Texture3D scattering_density_texture,
	SamplerState atmo_sampler,
	float r, float mu, float mu_s, float nu,
	bool ray_r_mu_intersects_ground)
{
    const int SAMPLE_COUNT = 50;
    float dx = DistanceToNearestAtmosphereBoundary(atmosphere, r, mu, ray_r_mu_intersects_ground) / float(SAMPLE_COUNT);
    float3 rayleigh_mie_sum = (float3) (0.0);
    for (int i = 0; i <= SAMPLE_COUNT; ++i)
    {
        float d_i = float(i) * dx;
        float r_i = ClampRadius(atmosphere, sqrt(d_i * d_i + 2.0 * r * mu * d_i + r * r));
        float mu_i = ClampCosine((r * mu + d_i) / r_i);
        float mu_s_i = ClampCosine((r * mu_s + d_i * nu) / r_i);
        float3 rayleigh_mie_i = GetScattering(atmosphere, scattering_density_texture, atmo_sampler, r_i, mu_i, mu_s_i, nu, ray_r_mu_intersects_ground) * GetTransmittance(atmosphere, transmittance_texture, atmo_sampler, r, mu, d_i, ray_r_mu_intersects_ground) * dx;
        float weight_i = (i == 0 || i == SAMPLE_COUNT) ? 0.5 : 1.0;
        rayleigh_mie_sum += rayleigh_mie_i * weight_i;
    }
    return rayleigh_mie_sum;
}

float3 ComputeScatteringDensityTexture(
	const in AtmosphereParameters atmosphere,
	Texture2D transmittance_texture,
	Texture3D single_rayleigh_scattering_texture,
	Texture3D single_mie_scattering_texture,
	Texture3D multiple_scattering_texture,
	Texture2D irradiance_texture,
	SamplerState atmo_sampler,
	float3 frag_coord,
	int scattering_order)
{
    float r;
    float mu;
    float mu_s;
    float nu;
    bool ray_r_mu_intersects_ground;
    GetRMuMuSNuFromScatteringTextureFragCoord(atmosphere, frag_coord, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    return ComputeScatteringDensity(atmosphere, transmittance_texture, single_rayleigh_scattering_texture, single_mie_scattering_texture, multiple_scattering_texture, irradiance_texture, atmo_sampler, r, mu, mu_s, nu, scattering_order);
}

float3 ComputeMultipleScatteringTexture(
	const in AtmosphereParameters atmosphere,
	Texture2D transmittance_texture,
	Texture3D scattering_density_texture,
	SamplerState atmo_sampler,
	float3 frag_coord, out float nu)
{
    float r;
    float mu;
    float mu_s;
    bool ray_r_mu_intersects_ground;
    GetRMuMuSNuFromScatteringTextureFragCoord(atmosphere, frag_coord, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    return ComputeMultipleScattering(atmosphere, transmittance_texture, scattering_density_texture, atmo_sampler, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
}

float3 ComputeDirectIrradiance(
	const in AtmosphereParameters atmosphere,
	Texture2D transmittance_texture,
	SamplerState atmo_sampler,
	float r, float mu_s)
{
    float alpha_s = atmosphere.sun_angular_radius;
    float average_cosine_factor = (mu_s < -alpha_s) ? 0.0 : ((mu_s > alpha_s) ? mu_s : (mu_s + alpha_s) * (mu_s + alpha_s) / (4.0 * alpha_s));
    return atmosphere.solar_irradiance * GetTransmittanceToTopAtmosphereBoundary(atmosphere, transmittance_texture, atmo_sampler, r, mu_s) * average_cosine_factor;
}

float3 ComputeIndirectIrradiance(
	const in AtmosphereParameters atmosphere,
	Texture3D single_rayleigh_scattering_texture,
	Texture3D single_mie_scattering_texture,
	Texture3D multiple_scattering_texture,
	SamplerState atmo_sampler,
	float r, float mu_s, int scattering_order)
{
    const int SAMPLE_COUNT = 32;
    const float dphi = PI / float(SAMPLE_COUNT);
    const float dtheta = PI / float(SAMPLE_COUNT);
    float3 result = (float3) (0.0);
    float3 omega_s = float3(sqrt(1.0 - mu_s * mu_s), 0.0, mu_s);
    for (int j = 0; j < SAMPLE_COUNT / 2; ++j)
    {
        float theta = (float(j) + 0.5) * dtheta;
        for (int i = 0; i < 2 * SAMPLE_COUNT; ++i)
        {
            float phi = (float(i) + 0.5) * dphi;
            float3 omega = float3(cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));
            float domega = dtheta * dphi * sin(theta);
            float nu = dot(omega, omega_s);
            result += GetScattering(atmosphere, single_rayleigh_scattering_texture, single_mie_scattering_texture, multiple_scattering_texture, atmo_sampler, r, omega.z, mu_s, nu, false /* ray_r_theta_intersects_ground */, scattering_order) * omega.z * domega;
        }
    }
    return result;
}

float3 ComputeDirectIrradianceTexture(
	const in AtmosphereParameters atmosphere,
	Texture2D transmittance_texture,
	SamplerState atmo_sampler,
	float2 frag_coord)
{
    float r;
    float mu_s;
    GetRMuSFromIrradianceTextureUv(atmosphere, frag_coord / IRRADIANCE_TEXTURE_SIZE, r, mu_s);
    return ComputeDirectIrradiance(atmosphere, transmittance_texture, atmo_sampler, r, mu_s);
}

float3 ComputeIndirectIrradianceTexture(
	const in AtmosphereParameters atmosphere,
	Texture3D single_rayleigh_scattering_texture,
	Texture3D single_mie_scattering_texture,
	Texture3D multiple_scattering_texture,
	SamplerState atmo_sampler,
	float2 frag_coord, int scattering_order)
{
    float r;
    float mu_s;
    GetRMuSFromIrradianceTextureUv(atmosphere, frag_coord / IRRADIANCE_TEXTURE_SIZE, r, mu_s);
    return ComputeIndirectIrradiance(atmosphere, single_rayleigh_scattering_texture, single_mie_scattering_texture, multiple_scattering_texture, atmo_sampler, r, mu_s, scattering_order);
}

#ifdef COMBINED_SCATTERING_TEXTURES
float3 GetExtrapolatedSingleMieScattering(const in AtmosphereParameters atmosphere, float4 scattering)
{
    if (scattering.r == 0.0)
    {
        return (float3) (0.0);
    }
    return scattering.rgb * scattering.a / scattering.r * (atmosphere.rayleigh_scattering.r / atmosphere.mie_scattering.r) * (atmosphere.mie_scattering / atmosphere.rayleigh_scattering);
}
#endif

float3 GetCombinedScattering(
	const in AtmosphereParameters atmosphere,
	Texture3D scattering_texture,
	Texture3D single_mie_scattering_texture,
	SamplerState atmo_sampler,
	float r, float mu, float mu_s, float nu,
	bool ray_r_mu_intersects_ground,
	out float3 single_mie_scattering)
{
    float4 uvwz = GetScatteringTextureUvwzFromRMuMuSNu(atmosphere, r, mu, mu_s, nu, ray_r_mu_intersects_ground);
    float tex_coord_x = uvwz.x * float(SCATTERING_TEXTURE_NU_SIZE - 1);
    float tex_x = floor(tex_coord_x);
    float _lerp = tex_coord_x - tex_x;
    float3 uvw0 = float3((tex_x + uvwz.y) / float(SCATTERING_TEXTURE_NU_SIZE), uvwz.z, uvwz.w);
    float3 uvw1 = float3((tex_x + 1.0 + uvwz.y) / float(SCATTERING_TEXTURE_NU_SIZE), uvwz.z, uvwz.w);
#ifdef COMBINED_SCATTERING_TEXTURES
    float4 combined_scattering =
		scattering_texture.Sample(atmo_sampler, uvw0) * (1.0 - _lerp) +
		scattering_texture.Sample(atmo_sampler, uvw1) * _lerp;
    float3 scattering = combined_scattering.xyz;
    single_mie_scattering = GetExtrapolatedSingleMieScattering(atmosphere, combined_scattering);
#else
	float3 scattering = (
		scattering_texture.Sample(atmo_sampler, uvw0) * (1.0 - _lerp) +
		scattering_texture.Sample(atmo_sampler, uvw1) * _lerp).xyz;
	single_mie_scattering = (
		single_mie_scattering_texture.Sample(atmo_sampler, uvw0) * (1.0 - _lerp) +
		single_mie_scattering_texture.Sample(atmo_sampler, uvw1) * _lerp).xyz;
#endif
    return scattering;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// RADIANCE API
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float3 GetSkyRadiance(
	const in AtmosphereParameters atmosphere,
	Texture2D transmittance_texture,
	Texture3D scattering_texture,
	Texture3D single_mie_scattering_texture,
	SamplerState atmo_sampler,
	float3 camera, float3 view_ray, float shadow_length,
	float3 sun_direction, out float3 transmittance)
{
    float r = length(camera);
    float rmu = dot(camera, view_ray);
    float distance_to_top_atmosphere_boundary = -rmu - sqrt(rmu * rmu - r * r + atmosphere.top_radius * atmosphere.top_radius);
    if (distance_to_top_atmosphere_boundary > 0.0)
    {
        camera = camera + view_ray * distance_to_top_atmosphere_boundary;
        r = atmosphere.top_radius;
        rmu += distance_to_top_atmosphere_boundary;
    }
    else if (r > atmosphere.top_radius)
    {
        transmittance = (float3) (1.0);
        return (float3) (0.0);
    }

    float mu = rmu / r;
    float mu_s = dot(camera, sun_direction) / r;
    float nu = dot(view_ray, sun_direction);
    bool ray_r_mu_intersects_ground = RayIntersectsGround(atmosphere, r, mu);
    transmittance = ray_r_mu_intersects_ground ? (float3) (0.0) : GetTransmittanceToTopAtmosphereBoundary(atmosphere, transmittance_texture, atmo_sampler, r, mu);
    float3 single_mie_scattering;
    float3 scattering;
    if (shadow_length == 0.0)
    {
        scattering = GetCombinedScattering(atmosphere, scattering_texture, single_mie_scattering_texture, atmo_sampler, r, mu, mu_s, nu, ray_r_mu_intersects_ground, single_mie_scattering);
    }
    else
    {
        float d = shadow_length;
        float r_p = ClampRadius(atmosphere, sqrt(d * d + 2.0 * r * mu * d + r * r));
        float mu_p = (r * mu + d) / r_p;
        float mu_s_p = (r * mu_s + d * nu) / r_p;
        scattering = GetCombinedScattering(atmosphere, scattering_texture, single_mie_scattering_texture, atmo_sampler, r_p, mu_p, mu_s_p, nu, ray_r_mu_intersects_ground, single_mie_scattering);
        float3 shadow_transmittance = GetTransmittance(atmosphere, transmittance_texture, atmo_sampler, r, mu, shadow_length, ray_r_mu_intersects_ground);
        scattering = scattering * shadow_transmittance;
        single_mie_scattering = single_mie_scattering * shadow_transmittance;
    }

    //return float3(mu, mu_s, nu);
    return scattering * RayleighPhaseFunction(nu) + single_mie_scattering * MiePhaseFunction(atmosphere.mie_phase_function_g, nu);
}

float3 GetSkyRadianceToPoint(
	const in AtmosphereParameters atmosphere,
	Texture2D transmittance_texture,
	Texture3D scattering_texture,
	Texture3D single_mie_scattering_texture,
	SamplerState atmo_sampler,
	float3 camera, float3 pos, float shadow_length,
	float3 sun_direction, out float3 transmittance)
{
    float3 view_ray = normalize(pos - camera);
    float r = length(camera);
    float rmu = dot(camera, view_ray);
    float distance_to_top_atmosphere_boundary = -rmu - sqrt(rmu * rmu - r * r + atmosphere.top_radius * atmosphere.top_radius);
    if (distance_to_top_atmosphere_boundary > 0.0)
    {
        camera = camera + view_ray * distance_to_top_atmosphere_boundary;
        r = atmosphere.top_radius;
        rmu += distance_to_top_atmosphere_boundary;
    }

    float mu = rmu / r;
    float mu_s = dot(camera, sun_direction) / r;
    float nu = dot(view_ray, sun_direction);
    float d = length(pos - camera);
    bool ray_r_mu_intersects_ground = RayIntersectsGround(atmosphere, r, mu);
    transmittance = GetTransmittance(atmosphere, transmittance_texture, atmo_sampler, r, mu, d, ray_r_mu_intersects_ground);
    float3 single_mie_scattering;
    float3 scattering = GetCombinedScattering(atmosphere, scattering_texture, single_mie_scattering_texture, atmo_sampler, r, mu, mu_s, nu, ray_r_mu_intersects_ground, single_mie_scattering);
    d = max(d - shadow_length, 0.0);
    float r_p = ClampRadius(atmosphere, sqrt(d * d + 2.0 * r * mu * d + r * r));
    float mu_p = (r * mu + d) / r_p;
    float mu_s_p = (r * mu_s + d * nu) / r_p;
    float3 single_mie_scattering_p;
    float3 scattering_p = GetCombinedScattering(atmosphere, scattering_texture, single_mie_scattering_texture, atmo_sampler, r_p, mu_p, mu_s_p, nu, ray_r_mu_intersects_ground, single_mie_scattering_p);
    float3 shadow_transmittance = transmittance;
    if (shadow_length > 0.0)
    {
        shadow_transmittance = GetTransmittance(atmosphere, transmittance_texture, atmo_sampler, r, mu, d, ray_r_mu_intersects_ground);
    }
    scattering = scattering - shadow_transmittance * scattering_p;
    single_mie_scattering = single_mie_scattering - shadow_transmittance * single_mie_scattering_p;
#ifdef COMBINED_SCATTERING_TEXTURES
    single_mie_scattering = GetExtrapolatedSingleMieScattering(atmosphere, float4(scattering, single_mie_scattering.r));
#endif
    single_mie_scattering = single_mie_scattering * smoothstep(float(0.0), float(0.01), mu_s);
    return scattering * RayleighPhaseFunction(nu) + single_mie_scattering * MiePhaseFunction(atmosphere.mie_phase_function_g, nu);
}

float3 GetSunAndSkyIrradiance(
	const in AtmosphereParameters atmosphere,
	Texture2D transmittance_texture,
	Texture2D irradiance_texture,
	SamplerState atmo_sampler,
	float3 pos, float3 normal, float3 sun_direction,
	out float3 sky_irradiance)
{
    float r = length(pos);
    float mu_s = dot(pos, sun_direction) / r;
    sky_irradiance = GetIrradiance(atmosphere, irradiance_texture, atmo_sampler, r, mu_s) * (1.0 + dot(normal, pos) / r) * 0.5;
    return atmosphere.solar_irradiance * GetTransmittanceToSun(atmosphere, transmittance_texture, atmo_sampler, r, mu_s) * max(dot(normal, sun_direction), 0.0);
}

float3 GetSolarRadiance(AtmosphereParameters atmosphere)
{
    return atmosphere.solar_irradiance / (PI * atmosphere.sun_angular_radius * atmosphere.sun_angular_radius);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// LUMINANCE API
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

float3 GetSkyLuminance(
    AtmosphereParameters atmosphere,
	Texture2D transmittance_texture,
	Texture3D scattering_texture,
	Texture3D single_mie_scattering_texture,
	SamplerState atmo_sampler,
	float3 camera, float3 view_ray, float shadow_length,
	float3 sun_direction, out float3 transmittance)
{
    return GetSkyRadiance(atmosphere, transmittance_texture, scattering_texture, single_mie_scattering_texture, atmo_sampler, camera, view_ray, shadow_length, sun_direction, transmittance) * SKY_SPECTRAL_RADIANCE_TO_LUMINANCE;
}

float3 GetSkyLuminanceToPoint(
    AtmosphereParameters atmosphere,
	Texture2D transmittance_texture,
	Texture3D scattering_texture,
    Texture3D single_mie_scattering_texture,
	SamplerState atmo_sampler,
	float3 camera, float3 pos, float shadow_length,
	float3 sun_direction, out float3 transmittance)
{
    return GetSkyRadianceToPoint(atmosphere, transmittance_texture, scattering_texture, single_mie_scattering_texture, atmo_sampler, camera, pos, shadow_length, sun_direction, transmittance) * SKY_SPECTRAL_RADIANCE_TO_LUMINANCE;
}

float3 GetSunAndSkyIlluminance(
    AtmosphereParameters atmosphere,
	Texture2D transmittance_texture,
	Texture2D irradiance_texture,
	SamplerState atmo_sampler,
	float3 p, float3 normal, float3 sun_direction,
	out float3 sky_irradiance)
{
    float3 sun_irradiance = GetSunAndSkyIrradiance(atmosphere, transmittance_texture, irradiance_texture, atmo_sampler, p, normal, sun_direction, sky_irradiance);
    sky_irradiance *= SKY_SPECTRAL_RADIANCE_TO_LUMINANCE;
    return sun_irradiance * SUN_SPECTRAL_RADIANCE_TO_LUMINANCE;
}

float3 GetSolarLuminance(AtmosphereParameters atmosphere)
{
    return atmosphere.solar_irradiance / (PI * atmosphere.sun_angular_radius * atmosphere.sun_angular_radius) * SUN_SPECTRAL_RADIANCE_TO_LUMINANCE;
}
