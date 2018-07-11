#include "common.h"

#include "window.h"
#include "renderer.h"
#include "gui.h"
#include "error.h"

extern LRESULT ImGui_ImplWin32_WndProcHandler(HWND hWnd, UINT msg, WPARAM wParam, LPARAM lParam);

namespace window {

	static int window_width;
	static int window_height;
	static const char *p_window_class_name = "Atmosphere Demo DX11 Window Class";
	static const char *p_window_name = "Atmosphere Demo DX11";
	HWND h_window = NULL;

	LRESULT CALLBACK window_proc(HWND h_window, UINT msg, WPARAM w_param, LPARAM l_param)
	{
		if(ImGui_ImplWin32_WndProcHandler(h_window, msg, w_param, l_param)) {
			return true;
		}

		switch(msg) {
			case WM_SIZE: {
				if(renderer::resize(l_param)) return 0;
			} break;

			case WM_DESTROY: {
				PostQuitMessage(0);
			} break;

			default: return DefWindowProc(h_window, msg, w_param, l_param);
		}

		return 0;
	}

	void init(HINSTANCE hInstance, uint32_t width, uint32_t height) {
		window_width = width;
		window_height = height;

		WNDCLASSEX window_class = { 0 };
		window_class.cbSize = sizeof(WNDCLASSEX);
		window_class.style = CS_HREDRAW | CS_VREDRAW;
		window_class.lpfnWndProc = window_proc;
		window_class.hInstance = hInstance;
		window_class.hIcon = LoadIcon(NULL, IDI_APPLICATION);
		window_class.hCursor = LoadCursor(NULL, IDC_ARROW);
		window_class.hbrBackground = NULL;
		window_class.lpszClassName = p_window_class_name;

		ATOM result = RegisterClassExA(&window_class);
		if(!result) { error_win32("RegisterClassExA", GetLastError()); return; };

		RECT window_rect = { 0, 0, window_width, window_height };
		AdjustWindowRect(&window_rect, WS_OVERLAPPEDWINDOW, FALSE);

		h_window = CreateWindowExA(
			0, p_window_class_name, p_window_name,
			WS_OVERLAPPEDWINDOW | WS_SYSMENU, CW_USEDEFAULT, CW_USEDEFAULT,
			window_rect.right - window_rect.left,
			window_rect.bottom - window_rect.top,
			NULL, NULL, hInstance, NULL);
		if(!h_window) { error_win32("CreateWindowExA", GetLastError()); return; };

		ShowWindow(h_window, SW_SHOW);
		UpdateWindow(h_window);
	}

	void get_client_size(uint32_t &width, uint32_t &height){
		width = window_width;
		height = window_height;
	}

	HWND get_handle() {
		return h_window;
	}

} // namespace window