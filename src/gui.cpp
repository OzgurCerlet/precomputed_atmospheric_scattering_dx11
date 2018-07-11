#include "common.h"

#include "external/imgui.h"
#include "external/imgui_impl_dx11.h"

#include "gui.h"
#include "window.h"
#include "renderer.h"
#include "atmosphere.h"

namespace gui  {

	static GuiData gui_data;
	static ImVec2 last_mouse_pos;
	static const float mouse_pos_scale = 1.0f / 500.0f;

	void init() {
		ImGui::CreateContext();
		ImGui_ImplDX11_Init((void*)window::get_handle(), renderer::get_device().Get(), renderer::get_device_context().Get());
	}

	void update() {
		ImGui_ImplDX11_NewFrame();
		ImVec2 mouse_pos = ImGui::GetMousePos();
		ImGuiIO& io = ImGui::GetIO();
		bool is_mouse_captured = io.WantCaptureMouse;
		bool is_ctrl_pressed = io.KeyCtrl;
		bool is_left_mouse_button_pressed = ImGui::IsMouseDown(0);
		if(!is_mouse_captured && is_left_mouse_button_pressed) {
			if(is_ctrl_pressed) {
				float sun_zenith_angle_in_radians = XMConvertToRadians(gui_data.sun_zenith_angle_in_degrees);
				sun_zenith_angle_in_radians -= (last_mouse_pos.y - mouse_pos.y) * mouse_pos_scale;
				sun_zenith_angle_in_radians = XMMax(0.f, XMMin(XM_PI, sun_zenith_angle_in_radians));
				gui_data.sun_zenith_angle_in_degrees = XMConvertToDegrees(sun_zenith_angle_in_radians);
				float sun_azimuth_angle_in_radians = XMConvertToRadians(gui_data.sun_azimuth_angle_in_degrees);
				sun_azimuth_angle_in_radians += (last_mouse_pos.x - mouse_pos.x) * mouse_pos_scale;
				gui_data.sun_azimuth_angle_in_degrees = XMConvertToDegrees(sun_azimuth_angle_in_radians);
			}
			else {
				float view_zenith_angle_in_radians = XMConvertToRadians(gui_data.view_zenith_angle_in_degrees);
				view_zenith_angle_in_radians += (last_mouse_pos.y - mouse_pos.y) * mouse_pos_scale;
				view_zenith_angle_in_radians = XMMax(0.f, XMMin(XM_PIDIV2, view_zenith_angle_in_radians));
				gui_data.view_zenith_angle_in_degrees = XMConvertToDegrees(view_zenith_angle_in_radians);
				float view_azimuth_angle_in_radians = XMConvertToRadians(gui_data.view_azimuth_angle_in_degrees);
				view_azimuth_angle_in_radians += (last_mouse_pos.x - mouse_pos.x) * mouse_pos_scale;
				gui_data.view_azimuth_angle_in_degrees = XMConvertToDegrees(view_azimuth_angle_in_radians);
			}	
		}
		last_mouse_pos = mouse_pos;
		if(io.MouseWheel < 0.0) { gui_data.view_distance *= 1.05f; }
		else if(io.MouseWheel > 0.0) { gui_data.view_distance /= 1.05f; };

		// Prepare main window
		{	
			ImVec2 window_pos = ImVec2(io.DisplaySize.x - 100, 0);
			ImGui::SetNextWindowPos(ImVec2(0,0), ImGuiCond_Always, ImVec2(0, 0));
			ImGui::SetNextWindowBgAlpha(0.3f); // Transparent background

			ImGui::Begin("Controls", NULL, ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoSavedSettings);
			ImGui::PushItemWidth(100);
			ImGui::Text("View: ");
			ImGui::InputFloat("View azimuth angle", &gui_data.view_azimuth_angle_in_degrees);   
			ImGui::InputFloat("View zenith angle", &gui_data.view_zenith_angle_in_degrees);
			ImGui::InputFloat("Sun azimuth angle", &gui_data.sun_azimuth_angle_in_degrees);
			ImGui::InputFloat("Sun zenith angle", &gui_data.sun_zenith_angle_in_degrees);
			ImGui::InputFloat("View distance", &gui_data.view_distance);
			ImGui::InputFloat("Exposure", &gui_data.exposure, 0.f, 0.f, "%.6f");
			ImGui::SliderInt("Predefined views", &gui_data.predefined_view_index, 0, 9);
			ImGui::Checkbox("Do white balance", &gui_data.do_white_balance);
			ImGui::Separator();
			ImGui::Text("Atmosphere Model: ");
			ImGui::Checkbox("Use ozone layer", &gui_data.use_ozone_layer);
			ImGui::Checkbox("Use half precision", &gui_data.use_half_precision);
			const char* solar_spectrum_options[] = { "CONSTANT", "REALISTIC"};
			static int solar_spectrum_option_current = gui_data.use_constant_solar_spectrum ? 0 : 1;
			ImGui::Combo("Solar spectrum", &solar_spectrum_option_current, solar_spectrum_options, IM_ARRAYSIZE(solar_spectrum_options));
			if(solar_spectrum_option_current) { gui_data.use_constant_solar_spectrum = false;}
			else { gui_data.use_constant_solar_spectrum = true;}
			const char* luminance_options[] = { "OFF", "APPROXIMATE", "PRECOMPUTED" };
			ImGui::Combo("Use Luminance", (int*)&gui_data.use_luminance, luminance_options, IM_ARRAYSIZE(luminance_options));
			ImGui::Separator();
			//if(ImGui::Button("Debug Dump Textures")) { gui_data.debug_dump_textures = true; };
			ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / io.Framerate, io.Framerate);
			ImGui::End();			
		}	
	}

	void render_frame() {
		ImGui::Render();
		ImGui_ImplDX11_RenderDrawData(ImGui::GetDrawData());
	}

	void clean_up() {
		ImGui_ImplDX11_Shutdown();
		ImGui::DestroyContext();
	}

	GuiData& get_data() {
		return gui_data;
	}

} // namespace gui


