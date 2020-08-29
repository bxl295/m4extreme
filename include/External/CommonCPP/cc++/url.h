// Copyright (C) 2001 Open Source Telecom Corporation.
//  
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software 
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
// 
// As a special exception to the GNU General Public License, permission is 
// granted for additional uses of the text contained in its release 
// of Common C++.
// 
// The exception is that, if you link the Common C++ library with other files
// to produce an executable, this does not by itself cause the
// resulting executable to be covered by the GNU General Public License.
// Your use of that executable is in no way restricted on account of
// linking the Common C++ library code into it.
//
// This exception does not however invalidate any other reasons why
// the executable file might be covered by the GNU General Public License.
//
// This exception applies only to the code released under the 
// name Common C++.  If you copy code from other releases into a copy of
// Common C++, as the General Public License permits, the exception does
// not apply to the code that you add in this way.  To avoid misleading
// anyone as to the status of such modified files, you must delete
// this exception notice from them.
// 
// If you write modifications of your own for Common C++, it is your choice
// whether to permit this exception to apply to your modifications.
// If you do not wish that, delete this exception notice.  

#ifndef	__CCXX_URL_H__
#define	__CCXX_URL_H__

#ifndef CCXX_CONFIG_H_
#include <cc++/config.h>
#endif

#ifndef __CCXX_SOCKET_H__
#include <cc++/socket.h>
#endif

#ifdef	CCXX_NAMESPACES
namespace ost {
#endif

typedef enum
{
	URL_SUCCESS = 0,
	URL_UNREACHABLE,
	URL_MISSING,
	URL_DENIED,
	URL_INVALID,
	URL_FORBIDDEN,
	URL_UNAUTHORIZED,
	URL_RELOCATED,
	URL_FAILURE,
	URL_TIMEOUT
} urlerror_t;

typedef enum
{
	URL_ANONYMOUS_AUTH = 0,
	URL_BASIC_AUTH
} urlauth_t;

typedef enum
{
	URL_BINARY_ENCODING = 0,
	URL_CHUNKED_ENCODING
} urlencoding_t;

typedef	enum
{
	HTTP_GET_METHOD,
	HTTP_PUT_METHOD,
	HTTP_POST_METHOD,
	FTP_GET_METHOD,
	FTP_PUT_METHOD,
	FILE_GET_METHOD,
	FILE_PUT_METHOD
} urlmethod_t;

typedef	enum
{
	HTTP_PROTOCOL_1_0,
	HTTP_PROTOCOL_1_1
} urlprotocol_t;

class CCXX_CLASS_EXPORT URLStream;

/**
 * A URL processing version of TCPStream.
 *
 * @author David Sugar <dyfet@ostel.com>
 * @short C++ url processing stream class.
 */
class CCXX_CLASS_EXPORT URLStream : public TCPStream
{
private:
	const char *agent, *referer, *cookie, *pragma, *user, *password;
	InetHostAddress proxyHost, prevHost;
	tpport_t proxyPort, prevPort;
	InetHostAddress refhost;
	urlmethod_t urlmethod;
	urlencoding_t encoding;
	urlprotocol_t protocol;
	urlauth_t auth;
	timeout_t timeout;
	bool persistent;
	bool follow;
	unsigned chunk;
	char authbuf[128];

	urlerror_t getHTTPHeaders(char *buffer, unsigned bufsize);

protected:
	urlerror_t sendHTTPHeader(const char *path, const char **vars, int bufsize);
	int underflow(void);

	/**
	 * Derived method to receive and parse http "headers".
	 *
	 * @param header keyword.
	 * @param header keyword value.
	 */
	virtual void httpHeader(const char *header, const char *value)
		{return;};

public:
	/**
	 * Construct an instance of URL stream.
	 */
	URLStream();

	/**
	 * Line parsing with conversion.
	 *
	 * @return URLStream object.
	 * @param buffer to store.
	 * @param maximum buffer size.
	 */
	URLStream &getline(char *buffer, size_t len);

	/**
	 * Get URL data from a named stream of a known buffer size.
	 *
	 * @return url error code.
	 * @param name of resource.
	 * @param size of buffer.
	 */
	urlerror_t Get(const char *path, int buffer = 512);

	/**
	 * Submit URL with vars passed as argument array.  This submit
	 * assumes "GET" method.  Use "Post" member to perform post.
	 *
	 * @return url error code.
	 * @param name of resource.
	 * @param vars to set.
	 * @param size of buffer.
	 */
	urlerror_t Submit(const char *path, const char **vars, int buffer = 512);

	/**
	 * Post URL vars with post method.
	 *
	 * @return url error code.
	 * @param name of resource.
	 * @param vars to set.
	 * @param size of buffer.
	 */
	urlerror_t Post(const char *path, const char **vars, int buffer = 512);

	/**
	 * Used to fetch header information for a resource.
	 *
	 * @return url error code.
	 * @param anme of resource.
	 * @param size of buffer.
	 */
	urlerror_t Head(const char *path, int buffer = 512);

	/**
	 * Close the URL stream for a new connection.
	 */
	void Close();

	/**
	 * Set the referer url.
	 *
	 * @param referer string.
	 */
	void setReferer(const char *str);

	/**
	 * Set the cookie to pass.
	 *
	 * @param cookie string.
	 */
	inline void setCookie(const char *str)
		{cookie = str;};

	/**
	 * Set user id for the url.
	 *
	 * @param user id.
	 */
	inline void setUser(const char *str)
		{user = str;};

	/**
	 * Set password for the url.
	 *
	 * @param password.
	 */
	inline void setPassword(const char *str)
		{password = str;};

	/**
	 * Set authentication type for the url.
	 *
	 * @param authentication.
	 */
	void setAuthentication(urlauth_t a, const char *str = NULL);

	/**
	 * Set the pragmas.
	 *
	 * @param pragma setting.
	 */
	inline void setPragma(const char *str)
		{pragma = str;};

	/**
	 * Set the proxy server used.
	 *
	 * @param proxy host.
	 * @param proxy port.
	 */
	void setProxy(const char *host, tpport_t port);

	/**
	 * Set the agent.
	 *
	 * @param agent value.
	 */
	inline void setAgent(const char *str)
		{agent = str;};

	/**
	 * Get url method (and protocol) employed.
	 *
	 * @return url method in effect.
	 */
	inline urlmethod_t getMethod(void)
		{return urlmethod;};

	/**
	 * Set socket timeout characteristics for processing URL
	 * requests.  Set to 0 for no default timeouts.
	 *
	 * @param timeout to set.
	 */
	inline void setTimeout(timeout_t to)
		{timeout = to;};

	/**
	 * Specify url following.  Set to false to disable following
	 * of relocation requests.
	 *
	 * @param true to enable following.
	 */
	inline void setFollow(bool enable)
		{follow = enable;};

	/**
	 * Specify http protocol level being used.
	 *
	 * @param protocol level.
	 */
	inline void setProtocol(urlprotocol_t pro)
		{protocol = pro;};
};

extern CCXX_EXPORT(char*) urldecode(char *source, char *dest = NULL);
extern CCXX_EXPORT(char*) urlencode(const char *source, char *dest, unsigned size);

#ifdef	CCXX_NAMESPACES
};
#endif

#endif
/** EMACS **
 * Local variables:
 * mode: c++
 * c-basic-offset: 8
 * End:
 */
