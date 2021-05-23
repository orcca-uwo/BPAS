
#ifndef _TASK_SCHEDULER_HPP_
#define _TASK_SCHEDULER_HPP_

#ifndef TASK_SCHED_PARALLEL
	#if defined(SERIAL) && SERIAL
	#define TASK_SCHED_PARALLEL 0
	#else
	#define TASK_SCHED_PARALLEL 1
	#endif
#endif

#include <functional>

#if TASK_SCHED_PARALLEL
#include <thread>
#include "ExecutorThreadPool.hpp"
#else
#include <deque>
#endif

/**
 * A class to encapsulate sharing and scheduling tasks across
 * many processors.
 */
class TaskScheduler {

#if TASK_SCHED_PARALLEL
	ExecutorThreadPool& m_pool;
#else
	std::deque<std::function<void()>> tasks;
#endif


public:

	/**
	 * Construct a new TaskScheduler with a (potentially) underlying thread pool.
	 * That thread pool may be shared with AyncGeneratorPool objects
	 * so that tasks scheduled through this scheduler have higher priority.
	 *
	 */
#if TASK_SCHED_PARALLEL
	TaskScheduler() : m_pool(ExecutorThreadPool::getThreadPool()) {}
#else
	TaskScheduler() : tasks() {}
#endif

	~TaskScheduler() {}

	/**
	 * Add a new high priority task (function to execute) to be scheduled for execution.
	 * @param f: the function to execute
	 * @params args: the arguments to pass to the function on execution.
	 */
	template <class Function, class... Args>
	void addPriorityTask(Function&& f, Args... args) {
		std::function<void()> boundF = std::bind(std::forward<Function>(f), std::forward<Args>(args)...);
#if TASK_SCHED_PARALLEL
		m_pool.addPriorityTask(boundF);
#else
		tasks.push_front(std::move(boundF));
#endif
	}

	/**
	 * Add a new task (function to execute) to be scheduled for execution.
	 * @param f: the function to execute
	 * @params args: the arguments to pass to the function on execution.
	 */
	template <class Function, class... Args>
	void addTask(Function&& f, Args... args) {
		std::function<void()> boundF = std::bind(std::forward<Function>(f), std::forward<Args>(args)...);
#if TASK_SCHED_PARALLEL
		m_pool.addTask(boundF);
#else
		tasks.push_back(std::move(boundF));
#endif
	}

	/**
	 * A waiting function to block the caller's execution until all scheduled tasks have completed.
	 */
	void waitForAllTasks() {
#if TASK_SCHED_PARALLEL
		m_pool.waitForAllThreads();
#else
		while (!tasks.empty()) {
			std::function<void()> f = tasks.front();
			tasks.pop_front();
			f();
		}
#endif
	}

};


#endif
