#pragma once

void error(const char *p_func_name, const char *p_message);

void error(const char *p_func_name);

void error_win32(const char* p_func_name, DWORD last_error);

void error_win32(const char* p_func_name, const Microsoft::WRL::ComPtr<ID3DBlob> &com_error_blob);