#include "common.h"

#include "gui.h"
#include "window.h"
#include "renderer.h"
#include "atmosphere.h"

#pragma comment(lib, "dxgi.lib")
#pragma comment(lib, "d3d11.lib")
#pragma comment(lib, "dxguid.lib")
#pragma comment(lib, "d3dcompiler.lib")

void update() {
	gui::update();
	atmosphere::update();
	renderer::update();
}

void render() {
	renderer::render_frame();
	gui::render_frame();
	renderer::present_frame();
}

void clean_up() {
	gui::clean_up();
}

int WINAPI wWinMain(
	HINSTANCE hInstance, 
	HINSTANCE prevInstance, 
	LPWSTR cmdLine, 
	int cmdShow)
{
	window::init(hInstance, 960, 540);
	renderer::init();
	gui::init();
	atmosphere::init();

	MSG msg = { 0 };
	while(msg.message != WM_QUIT) {
		if(PeekMessageA(&msg, 0, 0, 0, PM_REMOVE)) {
			TranslateMessage(&msg);
			DispatchMessageA(&msg);
		}
		else {
			update();
			render();
		}
	}
	clean_up();
	return 0;
}