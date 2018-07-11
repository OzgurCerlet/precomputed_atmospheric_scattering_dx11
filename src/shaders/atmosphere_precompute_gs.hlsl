
struct GsInput
{
    float4 pos_cs : SV_POSITION;
    float2 uv : TEXCOORD;
    uint slice_index : SLICE;
};

struct GsOutput
{
    float4 pos_cs : SV_POSITION;
    float2 uv : TEXCOORD;
    uint slice_index : SV_RENDERTARGETARRAYINDEX;
};

[maxvertexcount(3)]
void gs_main(triangle GsInput input[3], inout TriangleStream<GsOutput> output)
{
    [unroll]
	for (uint i = 0; i < 3; i++){
        GsOutput tri;
		tri.pos_cs = input[i].pos_cs;
        tri.uv = input[i].uv;
        tri.slice_index = input[i].slice_index;
        output.Append(tri);
    }
}