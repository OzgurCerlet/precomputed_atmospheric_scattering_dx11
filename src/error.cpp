#include "common.h"

#include <winerror.h>
#include "error.h"

void error(const char *p_func_name, const char *p_message) {
	char display_msg[512] = { 0 };

	strcat_s(display_msg, p_func_name);
	strcat_s(display_msg, " failed with error: ");
	strcat_s(display_msg, p_message);

	int choice = MessageBoxA(NULL, display_msg, NULL, MB_ABORTRETRYIGNORE | MB_ICONERROR);
	switch(choice) {
		case IDABORT: exit(-1); break;
		case IDRETRY: DebugBreak(); break;
		case IDIGNORE: return;
		default: return;
	}
}

void error(const char *p_func_name) {
	error(p_func_name, "UNKNOWN!?");
}

void error_win32(const char* p_func_name, DWORD last_error) {
	void* p_error_msg = nullptr;
	FormatMessageA(
		FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
		NULL,
		last_error,
		MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
		(char*)&p_error_msg,
		0,
		NULL);

	error(p_func_name, (char*)p_error_msg);
	LocalFree(p_error_msg);
}

void error_win32(const char* p_func_name, const ComPtr<ID3DBlob> &com_error_blob) {
	if(com_error_blob) {
		error(p_func_name, (const char*)(com_error_blob->GetBufferPointer()));
	}
	else {
		error_win32(p_func_name, GetLastError());
	}
}

