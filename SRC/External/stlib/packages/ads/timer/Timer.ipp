// -*- C++ -*-

#if !defined(__ads_timer_Timer_ipp__)
#error This file is an implementation detail of Timer.
#endif

namespace ads {

#ifdef __APPLE__

// Use the mach absolute time units. Information is available at
// http://developer.apple.com/mac/library/qa/qa2004/qa1398.html

inline
Timer::
Timer() :
   // Initialize with an invalid start time.
   _start(std::numeric_limits<uint64_t>::max()),
   _elapsed(0) {
   mach_timebase_info(&_timeBaseInfo);
}

inline
void
Timer::
start() {
   _start = mach_absolute_time();
}

inline
void
Timer::
stop() {
#ifdef DEBUG_stlib
   // Check that the clock was started.
   assert(_start != std::numeric_limits<uint64_t>::max());
#endif
   _elapsed += mach_absolute_time() - _start;
}

inline
void
Timer::
reset() {
   _start = std::numeric_limits<uint64_t>::max();
   _elapsed = 0;
}

inline
Timer::
operator double() const {
   return 1e-9 *(_elapsed * _timeBaseInfo.numer / _timeBaseInfo.denom);
}

#else
//----------------------------------------------------------------------------

inline
Timer::
Timer() :
   // Initialize with an invalid start time.
   _start(std::numeric_limits<clock_t>::max()),
   _elapsed(0) {
}

inline
void
Timer::
start() {
   _start = std::clock();
}

inline
void
Timer::
stop() {
#ifdef DEBUG_stlib
   // Check that the clock was started.
   assert(_start != std::numeric_limits<clock_t>::max());
#endif
   _elapsed += std::clock() - _start;
}

inline
void
Timer::
reset() {
   _start = std::numeric_limits<clock_t>::max();
   _elapsed = 0;
}

inline
Timer::
operator double() const {
   return double(_elapsed) / CLOCKS_PER_SEC;
}

#endif

} // namespace ads
