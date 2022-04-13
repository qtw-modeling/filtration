#ifndef __LOGGER_H__
#define __LOGGER_H__

#include <ostream>

#include <set>
#include <vector>
#include <fstream>

namespace logger {

enum LoggerLevel {
	DEBUG = 0,
	INFO = 1,
	WARNING = 2,
	ERROR = 3,
	FATAL = 4,
	LEVEL_MAX = 5,
};

class BroadcastStream;

class Logger {
	static Logger *_instance;
	std::vector<BroadcastStream *> levels;
	Logger();
	~Logger();
public:
	static Logger &getInstance();
	static void cleanup();
	BroadcastStream &getStream(LoggerLevel lvl);
	static std::string prefix(LoggerLevel lvl);
};

class BroadcastStream {
	std::set<std::ostream *> streams;
	BroadcastStream(const BroadcastStream &);
public:
	BroadcastStream();
	explicit BroadcastStream(std::ostream &p);
	template <typename T>
	BroadcastStream &operator <<(const T &x) {
		for (std::set<std::ostream *>::iterator i = streams.begin(); i != streams.end(); i++)
			(**i) << x;
		return *this;
	}
	BroadcastStream &operator <<(std::ostream &(*pf)(std::ostream &)) {
		for (std::set<std::ostream *>::iterator i = streams.begin(); i != streams.end(); i++)
			(**i) << pf;
		return *this;
	}
	void addStream(std::ostream &p);
	void removeStream(std::ostream &p);
};

enum SilentLogger { SILENT_LOGGER = 0 };

template<typename T>
static inline SilentLogger operator <<(SilentLogger foo, T) {
	return foo;
}

static inline SilentLogger operator <<(SilentLogger foo, std::ostream &(*)(std::ostream &)) {
	return foo;
}

#define ENABLE_LOGGING 1

#if ENABLE_LOGGING
# define LOG(level) (logger::Logger::getInstance().getStream(logger::level) << logger::Logger::prefix(logger::level))
#else
# define LOG(level)	logger::SILENT_LOGGER
#endif

}; /* namespace logger */

#endif
