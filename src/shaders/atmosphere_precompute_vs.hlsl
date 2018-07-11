
struct VsOutput
{
    float4 pos_cs : SV_POSITION;
    noperspective float2 uv : TEXCOORD;
    nointerpolation uint slice_index : SLICE;
};

VsOutput vs_main(uint vertex_id : SV_VERTEXID, uint inst_id : SV_INSTANCEID)
{
    float2 uv = float2((vertex_id << 1) & 2, vertex_id & 2);
    float4 pos_cs = float4(uv * float2(2.0f, -2.0f) + float2(-1.0f, 1.0f), 0.0f, 1.0f);
 
    VsOutput result = (VsOutput) 0;
    result.uv = float2(uv.x, uv.y);
    result.pos_cs = pos_cs;
    result.slice_index = inst_id;

    return result;
}