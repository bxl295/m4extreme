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
// of APE.
// 
// The exception is that, if you link the APE library with other files
// to produce an executable, this does not by itself cause the
// resulting executable to be covered by the GNU General Public License.
// Your use of that executable is in no way restricted on account of
// linking the APE library code into it.
//
// This exception does not however invalidate any other reasons why
// the executable file might be covered by the GNU General Public License.
//
// This exception applies only to the code released under the 
// name APE.  If you copy code from other releases into a copy of
// APE, as the General Public License permits, the exception does
// not apply to the code that you add in this way.  To avoid misleading
// anyone as to the status of such modified files, you must delete
// this exception notice from them.
// 
// If you write modifications of your own for APE, it is your choice
// whether to permit this exception to apply to your modifications.
// If you do not wish that, delete this exception notice.  

#ifndef	__CCXX_XML_H__
#define	__CCXX_XML_H__

#ifndef CCXX_CONFIG_H_
#include <cc++/config.h>
#endif

#ifndef	__CCXX_THREAD_H__
#include <cc++/thread.h>
#endif

#ifndef	__CCXX_SLOG_H__
#include <cc++/slog.h>
#endif

#ifdef	CCXX_NAMESPACES
namespace ost {
#endif

/**
 * This class impliments a basic XML stream parser that can be used to
 * examine an XML resource thru virtual I/O methods.  This class must
 * be derived into one that can impliment the physical I/O required to
 * parse actual data.  A mixer class using XMLStream and URLStream would
 * seem a likely combination for this purpose.
 *
 * @author David Sugar <dyfet@ostel.com>
 * @short XML Stream Parser (SAX)
 */
class XMLStream
{
public:
	/**
	 * May perform an open operation on behalf of a parsed resource.
	 * In some cases, the parser may be merged with a class that
	 * already has performed some kind of open, and this method can
	 * then be ignored.
	 *
	 * @return true if open is successful.
	 * @param resource passed to Parse methods.
	 */
	virtual bool Open(const char *resource)
		{return true;};

	/**
	 * May perform a close operation of an i/o source when the parser
	 * has completed operation.
	 */
	virtual void Close(void)
		{return;};

	/**
	 * Get error logging level.
	 *
	 * @return error logging level.
	 */
	virtual slog_level_t getLogging(void)
		{return SLOG_CRITICAL;};

	/**
	 * Virtual to receive embedded comments in an XML document being
	 * parsed.
	 *
	 * @param text comment extracted.
	 */
	virtual void comment(const unsigned char *text)
		{return;};

	/**
	 * Read method to aquire data for the parser.
	 *
	 * @return number of bytes actually read.
	 * @param buffer to read data into.	
	 * @param number of bytes to read.
	 */
	virtual int Read(unsigned char *buffer, int len) = 0;

	/**
	 * Virtual to receive character text extracted from the document
	 * in the current element.
	 *
	 * @param text received.
	 * @param length of text received.
	 */
	virtual void characters(const unsigned char *text, unsigned len) = 0;

	/**
	 * Identify start of document event.
	 */
	virtual void startDocument(void)
		{return;};

	/**
	 * Identify end of document event.
	 */
	virtual void endDocument(void)
		{return;};

	/**
	 * Identify start of an element in the document.
	 *
	 * @param name of element found.
	 * @param list of attributes extracted.
	 */
	virtual void startElement(const unsigned char *name, const unsigned char **attr) = 0;

	/**
	 * Identify end of an element in the document.
	 *
	 * @param name of element found.
	 */
	virtual void endElement(const unsigned char *name) = 0;

	/**
	 * Parse a resource as a stream thru the virtual read method.
	 *
	 * @return true if well formed document has been fully parsed.
	 * @param optional name of resource.
	 */
	bool Parse(const char *resource = NULL);
};

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
