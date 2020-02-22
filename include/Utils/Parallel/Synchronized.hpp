
#ifndef _SYNCHRONIZED_HPP_
#define _SYNCHRONIZED_HPP_


/**
 * A special little code block macro meant to emulate the
 * keyword "synchronized" in java. 
 * NOTE: you can *not* break or return in the middle of a synchronized block.
 *       it will cause a dealock where the lock is never unlocked.
 */
#define synchronized(m) for(std::unique_lock<std::recursive_mutex> lk(m); lk; lk.unlock())

/**
 * A special little code block macro meant to emulate the
 * keyword "synchronized" in java. 
 * NOTE: you can *not* break or return in the middle of a synchronized block.
 *       it will cause a dealock where the lock is never unlocked.
 */
#define synchronized_nonrecursive(m) for(std::unique_lock<std::mutex> lk(m); lk; lk.unlock())


#endif