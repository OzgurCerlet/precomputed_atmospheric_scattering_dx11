#pragma once

namespace window {
	
	void init(HINSTANCE hInstance, uint32_t width, uint32_t height);

	void get_client_size(uint32_t &width, uint32_t &height);

	HWND get_handle();

 } // namespace window