#pragma once

struct ID3D11Device;
struct ID3D11DeviceContext;
namespace atmosphere { enum Luminance; }

namespace gui {
	
	typedef struct GuiData {
		float view_zenith_angle_in_degrees = 84.2248f;
		float view_azimuth_angle_in_degrees = -5.729578f;
		float view_distance = 9000.0f;
		float sun_zenith_angle_in_degrees = 74.48451f;
		float sun_azimuth_angle_in_degrees = 166.1578f;
		float exposure = 10.f;
		int predefined_view_index = 0;
		bool do_white_balance = true;

		atmosphere::Luminance use_luminance;
		bool use_ozone_layer = true;
		bool use_constant_solar_spectrum = false;
		bool use_half_precision = false;

		bool debug_dump_textures = false;
	} GuiData;

	void init();
	
	void update();
	
	void render_frame();
	
	void clean_up();

	GuiData& get_data();

} // namespace gui