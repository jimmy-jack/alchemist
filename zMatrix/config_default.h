#ifndef _CONFIG_H
#define _CONFIG_H
#include <Windows.h>

//#define DEBUG

#if defined(DEBUG)
#define _InsightLog_(str) do{\
	char buf[1000];\
	sprintf_s(buf, "--------------------------------	[%s]<%s>(%d)%s	----------------------------------------\n", __FILE__, __func__, __LINE__, str);\
	OutputDebugString(buf);\
}while (0)

#define _log_(str)  cout << "[" << __FILE__<<"]" << "<" << __func__ << ">" << "(" << __LINE__ << ")" << str << endl;

#else
#define _log_(str)
#endif

#endif