#include "common.h"

#include <dxgi1_5.h>
#include <d3dcompiler.h>
#include <d3d11sdklayers.h>

//#include "D:/Octarine/OctarineImage/octarine_image.h"

#include "renderer.h"
#include "gui.h"
#include "window.h"
#include "atmosphere.h"
#include "error.h"

#define MAX_NUM_DXGI_ADAPTERS 8
#define _D3D_SET_DEBUG_NAME(A,B)	((A)->SetPrivateData(WKPDID_D3DDebugObjectName, sizeof((B)), (B)))
#define D3D_SET_DEBUG_NAME(A)		_D3D_SET_DEBUG_NAME(A,_STRINGIZE(A))

namespace renderer {

	ComPtr<ID3D11Device> com_device = nullptr;
	ComPtr<ID3D11DeviceContext> com_device_context = nullptr;
	ComPtr<IDXGISwapChain1> com_swap_chain = nullptr;
	ComPtr<ID3D11RenderTargetView> com_backbuffer_rtv = nullptr;
	ComPtr<ID3D11RasterizerState> com_rasterizer_state = nullptr;
	ComPtr<ID3D11SamplerState> com_sampler_state = nullptr;
	ComPtr<ID3D11DepthStencilState> com_depth_stencil_state = nullptr;
	ComPtr<ID3D11BlendState> com_blend_state_01 = nullptr;
	ComPtr<ID3D11BlendState> com_blend_state_0011 = nullptr;

	static uint32_t backbuffer_width;
	static uint32_t backbuffer_height;
	static bool is_full_precision_rgb_supported = false;
	static int last_predefined_view_index = 0;
	static float fov_y_angle_deg = 50.f;
	static float near_plane = 1.0;
	static UINT shader_compile_flags;
	
	struct AtmosphereConstantBuffer {
		struct {
			XMFLOAT4X4 view_from_clip;
			XMFLOAT4X4 world_from_view;
			XMFLOAT3 view_pos_ws;
			float sun_disk_size_x;
			XMFLOAT3 earth_center_pos_ws;
			float sun_disk_size_y;
			XMFLOAT3 sun_direction_ws;
			float    exposure;
			XMFLOAT3 white_point;
			int layer;
			XMFLOAT4X3 luminance_from_radiance;
			int scattering_order;
			float _pad[3];
		} data;
		ComPtr<ID3D11Buffer> com_cb = nullptr;
	} atmosphere_cb = {};

	void create_texture_2d(uint32_t width, uint32_t height, D3D11_SUBRESOURCE_DATA *p_init_data, DXGI_FORMAT format, Texture2D* p_out_texture) {
		p_out_texture->com_tex.Reset();
		p_out_texture->com_srv.Reset();
		p_out_texture->com_rtv.Reset();

		HRESULT h_result = 0;
		D3D11_TEXTURE2D_DESC desc = {};
		desc.Width = width;
		desc.Height = height;
		desc.MipLevels = 1;
		desc.ArraySize = 1;
		desc.Format = format;
		desc.SampleDesc.Count = 1;
		desc.Usage = D3D11_USAGE_DEFAULT;
		desc.BindFlags = D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_RENDER_TARGET;

		h_result = com_device->CreateTexture2D(&desc, p_init_data, p_out_texture->com_tex.GetAddressOf());
		if(h_result != S_OK) error_win32("CreateTexture2D", (DWORD)h_result);

		D3D11_SHADER_RESOURCE_VIEW_DESC srv_desc = {};
		srv_desc.Format = format;
		srv_desc.ViewDimension = D3D11_SRV_DIMENSION_TEXTURE2D;
		srv_desc.Texture2D.MipLevels = 1;

		h_result = com_device->CreateShaderResourceView(p_out_texture->com_tex.Get(), &srv_desc, p_out_texture->com_srv.GetAddressOf());
		if(h_result != S_OK) error_win32("CreateShaderResourceView", (DWORD)h_result);

		D3D11_RENDER_TARGET_VIEW_DESC rtv_desc = {};
		rtv_desc.Format = format;
		rtv_desc.ViewDimension = D3D11_RTV_DIMENSION_TEXTURE2D;

		h_result = com_device->CreateRenderTargetView(p_out_texture->com_tex.Get(), &rtv_desc, p_out_texture->com_rtv.GetAddressOf());
		if(h_result != S_OK) error_win32("CreateRenderTargetView", (DWORD)h_result);
	}

	void create_texture_3d(uint32_t width, uint32_t height, uint32_t depth, D3D11_SUBRESOURCE_DATA *p_init_data, DXGI_FORMAT format, Texture3D* p_out_texture) {
		p_out_texture->com_tex.Reset();
		p_out_texture->com_srv.Reset();
		p_out_texture->com_rtv.Reset();

		HRESULT h_result = 0;
		D3D11_TEXTURE3D_DESC desc = {};
		desc.Width = width;
		desc.Height = height;
		desc.Depth = depth;
		desc.MipLevels = 1;
		desc.Format = format;
		desc.Usage = D3D11_USAGE_DEFAULT;
		desc.BindFlags = D3D11_BIND_SHADER_RESOURCE | D3D11_BIND_RENDER_TARGET;

		h_result = com_device->CreateTexture3D(&desc, p_init_data, p_out_texture->com_tex.GetAddressOf());
		if(h_result != S_OK) error_win32("CreateTexture3D", (DWORD)h_result);

		D3D11_SHADER_RESOURCE_VIEW_DESC srv_desc = {};
		srv_desc.Format = format;
		srv_desc.ViewDimension = D3D11_SRV_DIMENSION_TEXTURE3D;
		srv_desc.Texture3D.MipLevels = 1;

		h_result = com_device->CreateShaderResourceView(p_out_texture->com_tex.Get(), &srv_desc, p_out_texture->com_srv.GetAddressOf());
		if(h_result != S_OK) error_win32("CreateShaderResourceView", (DWORD)h_result);

		D3D11_RENDER_TARGET_VIEW_DESC rtv_desc = {};
		rtv_desc.Format = format;
		rtv_desc.ViewDimension = D3D11_RTV_DIMENSION_TEXTURE3D;
		rtv_desc.Texture3D.WSize = -1;

		h_result = com_device->CreateRenderTargetView(p_out_texture->com_tex.Get(), &rtv_desc, p_out_texture->com_rtv.GetAddressOf());
		if(h_result != S_OK) error_win32("CreateRenderTargetView", (DWORD)h_result);
	}
/*
	void debug_dump_texture_2d(const Texture2D& tex, const char *p_name) {

		ComPtr<ID3D11Texture2D> com_staging_tex = nullptr;

		HRESULT h_result = 0;
		D3D11_TEXTURE2D_DESC desc = {};
		tex.com_tex->GetDesc(&desc);
		desc.Usage = D3D11_USAGE_STAGING;
		desc.BindFlags = 0;
		desc.CPUAccessFlags = D3D11_CPU_ACCESS_READ;

		h_result = com_device->CreateTexture2D(&desc, nullptr, com_staging_tex.GetAddressOf());
		if(h_result != S_OK) error_win32("CreateTexture2D", (DWORD)h_result);

		com_device_context->CopyResource(com_staging_tex.Get(), tex.com_tex.Get());
		com_device_context->Flush();

		D3D11_MAPPED_SUBRESOURCE mapped_resource = {};
		h_result = com_device_context->Map(com_staging_tex.Get(), 0, D3D11_MAP_READ, 0, &mapped_resource);
		if(h_result != S_OK) error_win32("Map", (DWORD)h_result);

		uint32_t texture_size = mapped_resource.RowPitch * desc.Height;
		void *p_out_data = malloc(texture_size);

		//char a_debug_name[64] = {};
		//uint32_t tex_debug_name_len = count_of(a_debug_name);
		//tex.com_tex->GetPrivateData(WKPDID_D3DDebugObjectName, &tex_debug_name_len, (void*)a_debug_name);

		memcpy(p_out_data, mapped_resource.pData, texture_size);

		OctarineImageHeader header = {};
		header.magic_num = OCTARINE_MAGIC_CONST;
		header.width = desc.Width;
		header.height = desc.Height;
		header.depth = 1;
		header.array_size = 1;
		header.mip_levels = 1;
		header.flags = 0;
		header.format = octarine_image_make_format_from_dxgi_format(desc.Format);
		header.size_of_data = texture_size;
		OCTARINE_IMAGE result = octarine_image_write_to_file(p_name, &header, p_out_data);
		if(result != OCTARINE_IMAGE_OK) { error("octarine_image_write_to_file"); };

		free(p_out_data);
		com_device_context->Unmap(com_staging_tex.Get(), 0);
	}

	void debug_dump_texture_3d(const Texture3D& tex, const char *p_name) {

		ComPtr<ID3D11Texture3D> com_staging_tex = nullptr;

		HRESULT h_result = 0;
		D3D11_TEXTURE3D_DESC desc = {};
		tex.com_tex->GetDesc(&desc);
		desc.Usage = D3D11_USAGE_STAGING;
		desc.BindFlags = 0;
		desc.CPUAccessFlags = D3D11_CPU_ACCESS_READ;

		h_result = com_device->CreateTexture3D(&desc, nullptr, com_staging_tex.GetAddressOf());
		if(h_result != S_OK) error_win32("CreateTexture3D", (DWORD)h_result);

		com_device_context->CopyResource(com_staging_tex.Get(), tex.com_tex.Get());
		com_device_context->Flush();

		D3D11_MAPPED_SUBRESOURCE mapped_resource = {};
		h_result = com_device_context->Map(com_staging_tex.Get(), 0, D3D11_MAP_READ, 0, &mapped_resource);
		if(h_result != S_OK) error_win32("Map", (DWORD)h_result);

		uint32_t texture_size = mapped_resource.RowPitch * desc.Height * desc.Depth;
		void *p_out_data = malloc(texture_size);

		//char a_debug_name[64] = {};
		//uint32_t tex_debug_name_len = count_of(a_debug_name);
		//tex.com_tex->GetPrivateData(WKPDID_D3DDebugObjectName, &tex_debug_name_len, (void*)a_debug_name);

		memcpy(p_out_data, mapped_resource.pData, texture_size);

		OctarineImageHeader header = {};
		header.magic_num = OCTARINE_MAGIC_CONST;
		header.width = desc.Width;
		header.height = desc.Height;
		header.depth = desc.Depth;
		header.array_size = 1;
		header.mip_levels = 1;
		header.flags = 0;
		header.format = octarine_image_make_format_from_dxgi_format(desc.Format);
		header.size_of_data = texture_size;
		OCTARINE_IMAGE result = octarine_image_write_to_file(p_name, &header, p_out_data);
		if(result != OCTARINE_IMAGE_OK) { error("octarine_image_write_to_file"); };

		free(p_out_data);
		com_device_context->Unmap(com_staging_tex.Get(), 0);
	}
	*/
	void compile_and_create_shader(ComPtr<ID3D11VertexShader>& com_vs, const wchar_t* p_name) {
		HRESULT h_result;
		ComPtr<ID3DBlob> com_shader_blob = nullptr;
		ComPtr<ID3DBlob> com_error_blob = nullptr;

		com_vs.Reset();
		h_result = D3DCompileFromFile(
			p_name,
			NULL,
			D3D_COMPILE_STANDARD_FILE_INCLUDE,
			"vs_main",
			"vs_5_0",
			shader_compile_flags,
			0,
			com_shader_blob.GetAddressOf(),
			com_error_blob.GetAddressOf()
		);
		if(h_result != S_OK) error_win32("D3DCompileFromFile", com_error_blob);

		h_result = com_device->CreateVertexShader(com_shader_blob->GetBufferPointer(), com_shader_blob->GetBufferSize(), NULL, com_vs.GetAddressOf());
		if(h_result != S_OK) error_win32("CreateVertexShader", (DWORD)h_result);
	}

	void compile_and_create_shader(ComPtr<ID3D11PixelShader>& com_ps, const wchar_t* p_name, const D3D_SHADER_MACRO *p_macros){
		HRESULT h_result;
		ComPtr<ID3DBlob> com_shader_blob = nullptr;
		ComPtr<ID3DBlob> com_error_blob = nullptr;

		com_ps.Reset();
		h_result = D3DCompileFromFile(
			p_name,
			p_macros,
			D3D_COMPILE_STANDARD_FILE_INCLUDE,
			"ps_main",
			"ps_5_0",
			shader_compile_flags,
			0,
			com_shader_blob.GetAddressOf(),
			com_error_blob.GetAddressOf()
		);
		if(h_result != S_OK) error_win32("D3DCompileFromFile", com_error_blob);

		h_result = com_device->CreatePixelShader(com_shader_blob->GetBufferPointer(), com_shader_blob->GetBufferSize(), NULL, com_ps.GetAddressOf());
		if(h_result != S_OK) error_win32("CreatePixelShader", (DWORD)h_result);
	}

	void compile_and_create_shader(ComPtr<ID3D11GeometryShader>& com_gs, const wchar_t* p_name) {
		HRESULT h_result;
		ComPtr<ID3DBlob> com_shader_blob = nullptr;
		ComPtr<ID3DBlob> com_error_blob = nullptr;

		com_gs.Reset();
		h_result = D3DCompileFromFile(
			p_name,
			NULL,
			D3D_COMPILE_STANDARD_FILE_INCLUDE,
			"gs_main",
			"gs_5_0",
			shader_compile_flags,
			0,
			com_shader_blob.GetAddressOf(),
			com_error_blob.GetAddressOf()
		);
		if(h_result != S_OK) error_win32("D3DCompileFromFile", com_error_blob);

		h_result = com_device->CreateGeometryShader(com_shader_blob->GetBufferPointer(), com_shader_blob->GetBufferSize(), NULL, com_gs.GetAddressOf());
		if(h_result != S_OK) error_win32("CreateGeometryShader", (DWORD)h_result);
	}

	bool check_full_precision_rgb_support() {
		return is_full_precision_rgb_supported;
	}

	void create_cb(ComPtr<ID3D11Buffer>& com_buffer, const void *p_data, uint32_t size) {
		D3D11_BUFFER_DESC cb_desc = {};
		cb_desc.ByteWidth = size;
		cb_desc.Usage = D3D11_USAGE_DYNAMIC;
		cb_desc.BindFlags = D3D11_BIND_CONSTANT_BUFFER;
		cb_desc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;

		D3D11_SUBRESOURCE_DATA cb_init_data = {};
		cb_init_data.pSysMem = p_data;
		HRESULT h_result = com_device->CreateBuffer(&cb_desc, &cb_init_data, com_buffer.GetAddressOf());
		if(h_result != S_OK) error_win32("CreateBuffer", (DWORD)h_result);
	}
	
	void update_cb(ComPtr<ID3D11Buffer>& com_buffer, const void *p_data, uint32_t data_size) {
		D3D11_MAPPED_SUBRESOURCE mapped_resource;
		HRESULT h_result = com_device_context->Map(com_buffer.Get(), 0, D3D11_MAP_WRITE_DISCARD, 0, &mapped_resource);
		if(h_result != S_OK) error_win32("Map", (DWORD)h_result);
		memcpy(mapped_resource.pData, p_data, data_size);
		com_device_context->Unmap(com_buffer.Get(), 0);
	}

	void init() {
		window::get_client_size(backbuffer_width, backbuffer_height);
#ifdef _DEBUG
		shader_compile_flags = D3DCOMPILE_DEBUG | D3DCOMPILE_SKIP_OPTIMIZATION | D3DCOMPILE_SKIP_VALIDATION;
#else
		// NOTE(cerlet): For some reason shader compilation is unusually slow! Trying to skip some work!
		shader_compile_flags = D3DCOMPILE_SKIP_OPTIMIZATION | D3DCOMPILE_SKIP_VALIDATION;
#endif
		
		{
			ComPtr<IDXGIFactory5> com_dxgi_factory = nullptr;
			ComPtr<IDXGIAdapter1> a_com_dxgi_adapters[MAX_NUM_DXGI_ADAPTERS] = { nullptr };
#ifdef  _DEBUG
			UINT flags = DXGI_CREATE_FACTORY_DEBUG;
#else
			UINT flags = 0;
#endif // _DEBUG
			HRESULT h_result = CreateDXGIFactory2(flags, IID_PPV_ARGS(&com_dxgi_factory));
			if(FAILED(h_result)) { error_win32("CreateDXGIFactory2", (DWORD)h_result); return; };

			UINT adapter_index = 0;
			for(;;) {
				h_result = com_dxgi_factory->EnumAdapters1(adapter_index, a_com_dxgi_adapters[adapter_index].GetAddressOf());
				if(FAILED(h_result)) {
					if((h_result == DXGI_ERROR_NOT_FOUND && adapter_index == 0) || h_result != DXGI_ERROR_NOT_FOUND) {
						error_win32("EnumAdapters1", (DWORD)h_result);
						return;
					}
					else break;
				}
				adapter_index++;
			}

			UINT device_creation_flags = D3D11_CREATE_DEVICE_SINGLETHREADED | D3D11_CREATE_DEVICE_DISABLE_GPU_TIMEOUT;
#ifdef  _DEBUG
			device_creation_flags |= D3D11_CREATE_DEVICE_DEBUG;
#endif // _DEBUG
			const D3D_FEATURE_LEVEL a_d3d_feauture_levels[] = { D3D_FEATURE_LEVEL_11_1 };

			h_result = D3D11CreateDevice(
				a_com_dxgi_adapters[0].Get(),
				D3D_DRIVER_TYPE::D3D_DRIVER_TYPE_UNKNOWN,
				NULL,
				device_creation_flags,
				a_d3d_feauture_levels,
				count_of(a_d3d_feauture_levels),
				D3D11_SDK_VERSION,
				com_device.GetAddressOf(),
				NULL,
				com_device_context.GetAddressOf()
			);
			//TEMP!
			for(uint32_t i = 0; i < adapter_index; ++i) {
				DXGI_ADAPTER_DESC desc = {};
				a_com_dxgi_adapters[i]->GetDesc(&desc);
				OutputDebugStringW(desc.Description);
				OutputDebugString("\n");
			}

			if(h_result != S_OK) { error_win32("D3D11CreateDevice", (DWORD)h_result); return; };

			UINT format_support_flags;
			DXGI_FORMAT full_precision_format = DXGI_FORMAT_R32G32B32_FLOAT;
			com_device->CheckFormatSupport(full_precision_format, &format_support_flags);
			is_full_precision_rgb_supported = ((format_support_flags & D3D11_FORMAT_SUPPORT_RENDER_TARGET) != 0);

			DXGI_SWAP_CHAIN_DESC1 swap_chain_desc = { 0 };
			swap_chain_desc.Width = backbuffer_width;
			swap_chain_desc.Height = backbuffer_height;
			swap_chain_desc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
			swap_chain_desc.Stereo = false;
			swap_chain_desc.SampleDesc.Count = 1;
			swap_chain_desc.SampleDesc.Quality = 0;
			swap_chain_desc.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
			swap_chain_desc.BufferCount = 2;
			swap_chain_desc.Scaling = DXGI_SCALING_STRETCH;
			swap_chain_desc.SwapEffect = DXGI_SWAP_EFFECT_FLIP_SEQUENTIAL;
			swap_chain_desc.AlphaMode = DXGI_ALPHA_MODE_IGNORE;
			swap_chain_desc.Flags = 0;

			h_result = com_dxgi_factory->CreateSwapChainForHwnd(com_device.Get(), window::get_handle(), &swap_chain_desc, NULL, NULL, com_swap_chain.GetAddressOf());
			if(h_result != S_OK) { error_win32("CreateSwapChainForHwnd", (DWORD)h_result); return; };

			ComPtr<ID3D11Texture2D> com_backbuffer_tex = nullptr;
			h_result = com_swap_chain->GetBuffer(0, IID_PPV_ARGS(com_backbuffer_tex.GetAddressOf()));
			if(h_result != S_OK) { error_win32("GetBuffer", (DWORD)h_result); return; };

			h_result = com_device->CreateRenderTargetView(com_backbuffer_tex.Get(), NULL, com_backbuffer_rtv.GetAddressOf());
			if(h_result != S_OK) { error_win32("CreateRenderTargetView", (DWORD)h_result); return; };
			D3D_SET_DEBUG_NAME(com_backbuffer_rtv);
		}

		{
			create_cb(atmosphere_cb.com_cb, &atmosphere_cb.data, sizeof(atmosphere_cb.data));

			{
				D3D11_RASTERIZER_DESC rasterizer_desc = {};
				rasterizer_desc.FillMode = D3D11_FILL_SOLID;
				rasterizer_desc.CullMode = D3D11_CULL_NONE;
				rasterizer_desc.FrontCounterClockwise = true;

				HRESULT h_result = com_device->CreateRasterizerState(&rasterizer_desc, com_rasterizer_state.GetAddressOf());
				if(h_result != S_OK) error_win32("CreateRasterizerState", (DWORD)h_result);
			}

			{
				D3D11_SAMPLER_DESC sampler_desc = {};
				sampler_desc.Filter = D3D11_FILTER_MIN_MAG_LINEAR_MIP_POINT;
				sampler_desc.AddressU = D3D11_TEXTURE_ADDRESS_CLAMP;
				sampler_desc.AddressV = D3D11_TEXTURE_ADDRESS_CLAMP;
				sampler_desc.AddressW = D3D11_TEXTURE_ADDRESS_CLAMP;

				HRESULT h_result = com_device->CreateSamplerState(&sampler_desc, com_sampler_state.GetAddressOf());
				if(h_result != S_OK) error_win32("CreateSamplerState", (DWORD)h_result);
			}

			{
				D3D11_DEPTH_STENCIL_DESC depth_stencil_desc = {};
				HRESULT h_result = com_device->CreateDepthStencilState(&depth_stencil_desc, com_depth_stencil_state.GetAddressOf());
				if(h_result != S_OK) error_win32("CreateRasterizerState", (DWORD)h_result);
			}
			
			{
				{
					D3D11_BLEND_DESC blend_state_desc = {};
					blend_state_desc.IndependentBlendEnable = true;
					blend_state_desc.RenderTarget[0].RenderTargetWriteMask = D3D11_COLOR_WRITE_ENABLE_ALL;
					blend_state_desc.RenderTarget[1].BlendEnable = true;
					blend_state_desc.RenderTarget[1].SrcBlend = D3D11_BLEND_ONE;
					blend_state_desc.RenderTarget[1].DestBlend = D3D11_BLEND_ONE;
					blend_state_desc.RenderTarget[1].BlendOp = D3D11_BLEND_OP_ADD;
					blend_state_desc.RenderTarget[1].SrcBlendAlpha = D3D11_BLEND_ONE;
					blend_state_desc.RenderTarget[1].DestBlendAlpha = D3D11_BLEND_ONE;
					blend_state_desc.RenderTarget[1].BlendOpAlpha = D3D11_BLEND_OP_ADD;
					blend_state_desc.RenderTarget[1].RenderTargetWriteMask = D3D11_COLOR_WRITE_ENABLE_ALL;
					HRESULT h_result = com_device->CreateBlendState(&blend_state_desc, com_blend_state_01.GetAddressOf());
				}
				{
					D3D11_BLEND_DESC blend_state_desc = {};
					blend_state_desc.IndependentBlendEnable = true;
					blend_state_desc.RenderTarget[0].RenderTargetWriteMask = D3D11_COLOR_WRITE_ENABLE_ALL;
					blend_state_desc.RenderTarget[1].RenderTargetWriteMask = D3D11_COLOR_WRITE_ENABLE_ALL;
					blend_state_desc.RenderTarget[2].BlendEnable = true;
					blend_state_desc.RenderTarget[2].SrcBlend = D3D11_BLEND_ONE;
					blend_state_desc.RenderTarget[2].DestBlend = D3D11_BLEND_ONE;
					blend_state_desc.RenderTarget[2].BlendOp = D3D11_BLEND_OP_ADD;
					blend_state_desc.RenderTarget[2].SrcBlendAlpha = D3D11_BLEND_ONE;
					blend_state_desc.RenderTarget[2].DestBlendAlpha = D3D11_BLEND_ONE;
					blend_state_desc.RenderTarget[2].BlendOpAlpha = D3D11_BLEND_OP_ADD;
					blend_state_desc.RenderTarget[2].RenderTargetWriteMask = D3D11_COLOR_WRITE_ENABLE_ALL;
					blend_state_desc.RenderTarget[3].BlendEnable = true;
					blend_state_desc.RenderTarget[3].SrcBlend = D3D11_BLEND_ONE;
					blend_state_desc.RenderTarget[3].DestBlend = D3D11_BLEND_ONE;
					blend_state_desc.RenderTarget[3].BlendOp = D3D11_BLEND_OP_ADD;
					blend_state_desc.RenderTarget[3].SrcBlendAlpha = D3D11_BLEND_ONE;
					blend_state_desc.RenderTarget[3].DestBlendAlpha = D3D11_BLEND_ONE;
					blend_state_desc.RenderTarget[3].BlendOpAlpha = D3D11_BLEND_OP_ADD;
					blend_state_desc.RenderTarget[3].RenderTargetWriteMask = D3D11_COLOR_WRITE_ENABLE_ALL;
					HRESULT h_result = com_device->CreateBlendState(&blend_state_desc, com_blend_state_0011.GetAddressOf());
				}
			}
		}
	}

	void set_view(gui::GuiData& gui_data, float view_distance, float view_zenith_angle_in_radians, float view_azimuth_angle_in_radians, float sun_zenith_angle_in_radians, float sun_azimuth_angle_in_radians,float exposure) {
		gui_data.view_distance = view_distance;
		gui_data.view_zenith_angle_in_degrees = XMConvertToDegrees(view_zenith_angle_in_radians);
		gui_data.view_azimuth_angle_in_degrees = XMConvertToDegrees(view_azimuth_angle_in_radians);
		gui_data.sun_zenith_angle_in_degrees = XMConvertToDegrees(sun_zenith_angle_in_radians);
		gui_data.sun_azimuth_angle_in_degrees = XMConvertToDegrees(sun_azimuth_angle_in_radians);
		gui_data.exposure = exposure;
	}

	void update() {
		gui::GuiData& gui_data = gui::get_data();
		if(last_predefined_view_index != gui_data.predefined_view_index) {
			float exposure_luminance_factor = gui_data.use_luminance != atmosphere::Luminance::NONE ? 1e-5f : 1.0f;
			switch(gui_data.predefined_view_index) {
				case 1: { set_view(gui_data, 9000.0f, 1.47f, 0.0f, 1.3f, 3.0f, 10.0f * exposure_luminance_factor); } break;
				case 2: { set_view(gui_data, 9000.0f, 1.47f, 0.0f, 1.564f, -3.0f, 10.0f * exposure_luminance_factor); } break;
				case 3: { set_view(gui_data, 7000.0f, 1.57f, 0.0f, 1.54f, -2.96f, 10.0f * exposure_luminance_factor); } break;
				case 4: { set_view(gui_data, 7000.0f, 1.57f, 0.0f, 1.328f, -3.044f, 10.0f * exposure_luminance_factor); } break;
				case 5: { set_view(gui_data, 9000.0f, 1.39f, 0.0f, 1.2f, 0.7f, 10.0f * exposure_luminance_factor); } break;
				case 6: { set_view(gui_data, 9000.0f, 1.5f, 0.0f, 1.628f, 1.05f, 200.0f * exposure_luminance_factor); } break;
				case 7: { set_view(gui_data, 7000.0f, 1.43f, 0.0f, 1.57f, 1.34f, 40.0f * exposure_luminance_factor); } break;
				case 8: { set_view(gui_data, 2.7e6f, 0.81f, 0.0f, 1.57f, 2.0f, 10.0f * exposure_luminance_factor); } break;
				case 9: { set_view(gui_data, 1.2e7f, 0.0f, 0.0f, 0.93f, -2.0f, 10.0f * exposure_luminance_factor); } break;
				default: break;
			}
			last_predefined_view_index = gui_data.predefined_view_index;
		}

		{
			float fov_y_angle_rad = XMConvertToRadians(fov_y_angle_deg);
			float aspect_ratio = (float)backbuffer_width / backbuffer_height;
			float scale_y = (float)(1.0 / tan(fov_y_angle_rad / 2.0));
			float scale_x = scale_y / aspect_ratio;
			XMFLOAT4X4 clip_from_view = { // left-handed reversed-z infinite projection
				scale_x, 0.0, 0.0, 0.0,
				0.0, scale_y, 0.0, 0.0,
				0.0, 0.0, 0.0, near_plane,
				0.0, 0.0, 1.0, 0.0
			};

			XMMATRIX view_from_clip = XMMatrixInverse(nullptr, XMLoadFloat4x4(&clip_from_view));
			XMStoreFloat4x4(&atmosphere_cb.data.view_from_clip, view_from_clip);

			float cos_theta = (float)cos(XMConvertToRadians(gui_data.view_zenith_angle_in_degrees));
			float sin_theta = (float)sin(XMConvertToRadians(gui_data.view_zenith_angle_in_degrees));
			float cos_phi = (float)cos(XMConvertToRadians(gui_data.view_azimuth_angle_in_degrees));
			float sin_phi = (float)sin(XMConvertToRadians(gui_data.view_azimuth_angle_in_degrees));

			XMFLOAT3 view_x_ws = { -sin_phi, cos_phi, 0.f };
			XMFLOAT3 view_y_ws = { -cos_theta * cos_phi, -cos_theta * sin_phi, sin_theta };
			XMFLOAT3 view_z_ws = { -sin_theta * cos_phi, -sin_theta * sin_phi, -cos_theta };
			XMFLOAT4X4 world_from_view = {
				view_x_ws.x, view_y_ws.x, view_z_ws.x, -view_z_ws.x * gui_data.view_distance / 1000.f,
				view_x_ws.y, view_y_ws.y, view_z_ws.y, -view_z_ws.y * gui_data.view_distance / 1000.f,
				view_x_ws.z, view_y_ws.z, view_z_ws.z, -view_z_ws.z * gui_data.view_distance / 1000.f,
				0.0, 0.0, 0.0, 1.0
			};
			atmosphere_cb.data.world_from_view = world_from_view;

			atmosphere_cb.data.view_pos_ws = { world_from_view(0,3),world_from_view(1,3),world_from_view(2,3) };
			atmosphere_cb.data.earth_center_pos_ws = { 0.f,0.f, -6360.0 };
			cos_theta = (float)cos(XMConvertToRadians(gui_data.sun_zenith_angle_in_degrees));
			sin_theta = (float)sin(XMConvertToRadians(gui_data.sun_zenith_angle_in_degrees));
			cos_phi = (float)cos(XMConvertToRadians(gui_data.sun_azimuth_angle_in_degrees));
			sin_phi = (float)sin(XMConvertToRadians(gui_data.sun_azimuth_angle_in_degrees));
			atmosphere_cb.data.sun_direction_ws = {
				cos_phi * sin_theta,
				sin_phi * sin_theta,
				cos_theta
			};

			// White Balance
			double white_point_r = 1.0;
			double white_point_g = 1.0;
			double white_point_b = 1.0;
			if(gui_data.do_white_balance) {
				atmosphere::compute_white_point(&white_point_r, &white_point_g, &white_point_b);
			}
			atmosphere_cb.data.white_point = { (float)white_point_r, (float)white_point_g , (float)white_point_b };

			const double sun_angular_radius_rad = 0.00935 / 2.0;
			atmosphere_cb.data.sun_disk_size_x = (float)tan(sun_angular_radius_rad);
			atmosphere_cb.data.sun_disk_size_y = (float)cos(sun_angular_radius_rad);
			atmosphere_cb.data.exposure = gui_data.exposure;
			update_cb(atmosphere_cb.com_cb, &atmosphere_cb.data, sizeof(atmosphere_cb.data));
		}

		//if(gui_data.debug_dump_textures) {
		//	debug_dump_texture_2d(transmittance_texture, "precomputed_final_transmittance");
		//	debug_dump_texture_2d(irradiance_texture, "precomputed_final_irradiance");
		//	debug_dump_texture_3d(scattering_texture, "precomputed_final_scattering");
		//	debug_dump_texture_3d(single_mie_scattering_texture, "precomputed_final_single_mie_scattering");
		//	gui_data.debug_dump_textures = false;
		//}
	}

	int resize(LPARAM lparam) {
		if(com_swap_chain) {
			backbuffer_width = (UINT)LOWORD(lparam);
			backbuffer_height = (UINT)HIWORD(lparam);
			HRESULT h_result;
			ComPtr<ID3D11Texture2D> com_backbuffer_tex = nullptr;
			h_result = com_swap_chain->GetBuffer(0, IID_PPV_ARGS(com_backbuffer_tex.GetAddressOf()));
			if(h_result != S_OK) { error_win32("GetBuffer", (DWORD)h_result); return 0; };
			com_backbuffer_tex.Reset();
			com_backbuffer_rtv.Reset();
			com_swap_chain->ResizeBuffers(0, backbuffer_width, backbuffer_height, DXGI_FORMAT_UNKNOWN, 0);
			h_result = com_swap_chain->GetBuffer(0, IID_PPV_ARGS(com_backbuffer_tex.GetAddressOf()));
			if(h_result != S_OK) { error_win32("GetBuffer", (DWORD)h_result); return 0; };
			h_result = com_device->CreateRenderTargetView(com_backbuffer_tex.Get(), NULL, com_backbuffer_rtv.GetAddressOf());
			if(h_result != S_OK) { error_win32("CreateRenderTargetView", (DWORD)h_result); return 0; };
			return 1;
		}
		return 0;
	}

	ComPtr<ID3D11Device>& get_device() {
		return com_device;
	}

	ComPtr<ID3D11DeviceContext>& get_device_context(){
		return com_device_context;
	}

	ComPtr<ID3D11SamplerState>& get_sampler(){
		return com_sampler_state;
	}

	ComPtr<ID3D11RasterizerState>& get_rasterizer_state(){
		return com_rasterizer_state;
	}

	ComPtr<ID3D11DepthStencilState>& get_depth_stencil_state(){
		return com_depth_stencil_state;
	}

	ComPtr<ID3D11BlendState>& get_blend_state_01(){
		return com_blend_state_01;
	}

	ComPtr<ID3D11BlendState>& get_blend_state_0011(){
		return com_blend_state_0011;
	}

	void render_frame()
	{
		float clear_color[4] = { 0.0,1.0,0.0,0.0 };
		com_device_context->ClearRenderTargetView(com_backbuffer_rtv.Get(), clear_color);

		com_device_context->IASetPrimitiveTopology(D3D11_PRIMITIVE_TOPOLOGY_TRIANGLESTRIP);

		com_device_context->VSSetConstantBuffers(0, 1, atmosphere_cb.com_cb.GetAddressOf());
		com_device_context->VSSetShader(atmosphere::get_demo_vs().Get(), 0, 0);

		D3D11_VIEWPORT viewport = { 0.0, 0.0, (float)backbuffer_width, (float)backbuffer_height, 0.0, 1.0 };
		com_device_context->RSSetViewports(1, &viewport);
		com_device_context->RSSetState(com_rasterizer_state.Get());

		com_device_context->PSSetConstantBuffers(0, 1, atmosphere_cb.com_cb.GetAddressOf());
		ID3D11SamplerState *a_ps_samplers[] = { com_sampler_state.Get() };
		com_device_context->PSSetSamplers(0, count_of(a_ps_samplers), a_ps_samplers);
		ID3D11ShaderResourceView *a_srvs[] = { 
			atmosphere::get_transmittance_texture().com_srv.Get(), 
			atmosphere::get_scattering_texture().com_srv.Get(), 
			atmosphere::get_single_mie_scattering_texture().com_srv.Get(), 
			atmosphere::get_irradiance_texture().com_srv.Get() };
		com_device_context->PSSetShaderResources(0, count_of(a_srvs), a_srvs);
		com_device_context->PSSetShader(atmosphere::get_demo_ps().Get(), 0, 0);

		com_device_context->OMSetDepthStencilState(com_depth_stencil_state.Get(), 0);
		com_device_context->OMSetRenderTargets(1, com_backbuffer_rtv.GetAddressOf(), nullptr);

		com_device_context->Draw(4, 0);

		ID3D11ShaderResourceView *const a_null_srvs[count_of(a_srvs)] = {};
		com_device_context->PSSetShaderResources(0, count_of(a_srvs), a_null_srvs);
	}

	void present_frame() {
		com_swap_chain->Present(1, 0);
	}

} // namespace renderer