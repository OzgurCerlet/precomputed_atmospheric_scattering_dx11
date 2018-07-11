#pragma pack_matrix(row_major)
cbuffer atmosphere_cb : register(b0)
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
    int layer;
    float3x3 luminance_from_radiance;
    int scattering_order;
};

struct VsOutput
{
    float4 pos_cs       : SV_POSITION;
    float2 uv           : TEXCOORD0;
    float3 view_ray_ws  : VIEW_RAY_WS;
};

VsOutput vs_main(uint vertex_id : SV_VERTEXID)
{
    float2 uv = float2(vertex_id % 2, vertex_id >> 1);
    float4 pos_cs = float4(uv * float2(2.0f, -2.0f) + float2(-1.0f, 1.0f), 0.0f, 1.0f);
    float4 view_ray_vs = mul(view_from_clip, pos_cs);
    float3 view_ray_ws = mul(world_from_view, float4(view_ray_vs.xyz, 0)).xyz;
  
    VsOutput result = (VsOutput)0;
    result.uv = uv;
    result.pos_cs = pos_cs;
    result.view_ray_ws = view_ray_ws;

    return result;
}