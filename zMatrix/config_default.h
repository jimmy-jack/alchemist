#ifndef _CONFIG_H
#define _CONFIG_H

#define DEBUG

#if defined(DEBUG)
#define _log_(str) do{\
cout << "[" << __FILE__<<"]" << "<" << __func__ << ">" << "(" << __LINE__ << ")" << str << endl;\
}while (0)
#define __log__(str) cout << "[" << __FILE__<<"]" << "<" << __func__ << ">" << "(" << __LINE__ << ")" << str << endl;
#else
#define _log_(str)
#endif

#endif