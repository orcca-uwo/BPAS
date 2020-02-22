
#ifndef _EVENT_THREAD_EXECEPTION_HPP
#define _EVENT_THREAD_EXECEPTION_HPP

/** 
 * Custom exception for handling exceptions on the EventTread.
 */
class EventThreadException: public std::runtime_error
{
public:
	EventThreadException(const char* c) : std::runtime_error(c) {}
	EventThreadException(const std::string& s) : std::runtime_error(s) {}
};


#endif
