#pragma once

namespace renderer {

	struct Texture2D
	{
		Microsoft::WRL::ComPtr<ID3D11Texture2D> com_tex = nullptr;
		Microsoft::WRL::ComPtr<ID3D11ShaderResourceView> com_srv = nullptr;
		Microsoft::WRL::ComPtr<ID3D11RenderTargetView> com_rtv = nullptr;
	};

	struct Texture3D
	{
		Microsoft::WRL::ComPtr<ID3D11Texture3D> com_tex = nullptr;
		Microsoft::WRL::ComPtr<ID3D11ShaderResourceView> com_srv = nullptr;
		Microsoft::WRL::ComPtr<ID3D11RenderTargetView> com_rtv = nullptr;
	};

	void init();

	void update();

	void render_frame();

	void present_frame();

	int resize(LPARAM lparam);

	Microsoft::WRL::ComPtr<ID3D11Device> &get_device();

	Microsoft::WRL::ComPtr<ID3D11DeviceContext> &get_device_context();

	Microsoft::WRL::ComPtr<ID3D11SamplerState> &get_sampler();

	Microsoft::WRL::ComPtr<ID3D11RasterizerState> &get_rasterizer_state();

	Microsoft::WRL::ComPtr<ID3D11DepthStencilState> &get_depth_stencil_state();

	Microsoft::WRL::ComPtr<ID3D11BlendState> &get_blend_state_01();

	Microsoft::WRL::ComPtr<ID3D11BlendState> &get_blend_state_0011();

	void create_texture_2d(uint32_t width, uint32_t height, D3D11_SUBRESOURCE_DATA *p_init_data, DXGI_FORMAT format, Texture2D* p_out_texture);
	
	void create_texture_3d(uint32_t width, uint32_t height, uint32_t depth, D3D11_SUBRESOURCE_DATA *p_init_data, DXGI_FORMAT format, Texture3D* p_out_texture);
	
	void create_cb(Microsoft::WRL::ComPtr<ID3D11Buffer>& com_buffer, const void *p_data, uint32_t size);

	void update_cb(Microsoft::WRL::ComPtr<ID3D11Buffer>& com_buffer, const void *p_data, uint32_t data_size);

	void compile_and_create_shader(Microsoft::WRL::ComPtr<ID3D11VertexShader>& com_vs, const wchar_t* p_name);

	void compile_and_create_shader(Microsoft::WRL::ComPtr<ID3D11PixelShader>& com_ps, const wchar_t* p_name, const D3D_SHADER_MACRO *p_macros = NULL);

	void compile_and_create_shader(Microsoft::WRL::ComPtr<ID3D11GeometryShader>& com_gs, const wchar_t* p_name);

	bool check_full_precision_rgb_support();

} // namespace renderer