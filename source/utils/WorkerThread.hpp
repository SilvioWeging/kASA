/***************************************************************************
*  Part of kASA: https://github.com/SilvioWeging/kASA
*  Inspired by https://github.com/log4cplus/ThreadPool/blob/master/ThreadPool.h
*
*  Copyright (C) 2020 Silvio Weging <silvio.weging@gmail.com>
*
*  Distributed under the Boost Software License, Version 1.0.
*  (See accompanying file LICENSE_1_0.txt or copy at
*  http://www.boost.org/LICENSE_1_0.txt)
**************************************************************************/

#pragma once

#include <vector>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>

using namespace std;

class WorkerThread {

	vector<thread> T; // is supposed to hold one thread
	int32_t ID = 0;
	vector<function<void(const int32_t&)>> tasks;
	condition_variable cv_worker;
	condition_variable cv_master;
	mutex taskMutex;
	bool stop = false, bStarted = false;


	////////////////////////////////////////////////////////
public:
	inline WorkerThread() {
		tasks.reserve(1);
		startWorker();
	}

	////////////////////////////////////////////////////////
	inline void setID(const uint32_t& id) {
		this->ID = id;
	}

	////////////////////////////////////////////////////////
	inline ~WorkerThread() {
		std::unique_lock<std::mutex> lock(taskMutex);
		stop = true;
		cv_worker.notify_one();
		cv_worker.wait(lock, [this] { return this->T.empty(); });
	}

	////////////////////////////////////////////////////////
	inline void pushTask(const function<void(const int32_t&)>& task) {
		tasks.push_back([task](const int32_t& id) { task(id); });
	}

	////////////////////////////////////////////////////////
	inline size_t printVecSize() {
		return tasks.size();
	}

	////////////////////////////////////////////////////////
	inline void startThread() {
		if (tasks.size()) {
			bStarted = true;
			cv_worker.notify_one();
		}
	}

	////////////////////////////////////////////////////////
	inline void setNumberOfTasks(const size_t& number) {
		tasks.reserve(number);
	}

	////////////////////////////////////////////////////////
	inline bool isEmpty() const {
		return tasks.empty();
	}

	////////////////////////////////////////////////////////
	inline void waitUntilFinished() {
		if (bStarted) {
			unique_lock<std::mutex> lock(this->taskMutex);
			this->cv_master.wait(lock,
				[this] { return this->tasks.empty(); });
		}
	}


private:
	////////////////////////////////////////////////////////
	inline void startWorker() {
		T.emplace_back(
			[this](){
				while (true) {
					{
						std::unique_lock<std::mutex> lock(this->taskMutex);
						this->cv_worker.wait(lock, [this]() {return !(this->tasks.empty()) || this->stop; });
					}

					if (this->stop) {
						tasks.clear();
						this->T[0].detach();
						this->T.pop_back();
						this->cv_worker.notify_one();
						return;
					}

					for (const auto& entry : tasks) {
						entry(this->ID);
					}

					tasks.clear();
					bStarted = false;
					std::unique_lock<std::mutex> lock(this->taskMutex);
					cv_master.notify_one();
				}
			}
		);
#if __GNUC__ && !defined(__llvm__)
		auto native_thread = T[0].native_handle();
		pthread_attr_t attr;
		pthread_getattr_np(native_thread, &attr);
		pthread_attr_setstacksize(&attr, 1024); // reduce stacksize
#endif
	}

};

#include <queue>

class WorkerQueue {
	vector<unique_ptr<thread> > pool; // this holds the threads
	priority_queue<int32_t, vector<int32_t>> availableCores;
	queue< pair<function<void(const int32_t&)>,int32_t> > tasks; // this holds the tasks
	condition_variable cv_worker;
	condition_variable cv_master;
	mutex taskMutex;
	bool stop = false, bStarted = false;

public:
	/////////////////////////////////////////////////////////
	inline WorkerQueue(const int32_t& iNumOfThreads) {
		for (int32_t i = 0; i < iNumOfThreads; ++i) {
			startWorker(i);
		}
	}

	////////////////////////////////////////////////////////
	inline ~WorkerQueue() {
		std::unique_lock<std::mutex> lock(taskMutex);
		stop = true;
		cv_worker.notify_all();
		cv_worker.wait(lock, [this] { 
			for (size_t i = 0; i < this->pool.size(); ++i) {
				if (this->pool[i] != nullptr) {
					return false;
				}
			}
			return true; 
		});
	}

	////////////////////////////////////////////////////////
	inline void pushTask(const function<void(const int32_t&)>& task, const int32_t& iNumOfThreads) {
		tasks.push(make_pair([task](const int32_t& ithreads) { task(ithreads); }, iNumOfThreads));
	}
	
	////////////////////////////////////////////////////////
	inline void start() {
		if (!tasks.empty()) {
			bStarted = true;
			cv_worker.notify_all();
		}
	}

	////////////////////////////////////////////////////////
	inline void waitUntilFinished() {
		if (bStarted) {
			unique_lock<std::mutex> lock(this->taskMutex);
			this->cv_master.wait(lock,
				[this] { return this->tasks.empty(); });
		}
		bStarted = false;
	}

private:
	////////////////////////////////////////////////////////
	inline void startWorker(const size_t& index) {
		pool.emplace_back( new thread(
			[this,index]() {
				while (true) {

					pair<function<void(const int32_t&)>,int32_t> task;
					int32_t iNumOfAvailThreads = 0;
					bool bQueueisEmpty = false;
					{
						std::unique_lock<std::mutex> lock(this->taskMutex);
						this->cv_worker.wait(lock, [this]() {return !(this->tasks.empty()) || this->stop; });

						if (this->stop) {
							this->pool[index]->detach();
							this->pool[index].reset();
							this->cv_worker.notify_one();
							return;
						}

						if (!this->tasks.empty()) {
							task = std::move(this->tasks.front());
							this->tasks.pop();
							bQueueisEmpty = this->tasks.empty();
							if (!availableCores.empty()) {
								iNumOfAvailThreads = this->availableCores.top();
								this->availableCores.pop();
							}
						}
						else {
							continue;
						}
					}

					if (bQueueisEmpty) {
						std::unique_lock<std::mutex> lock(this->taskMutex);
						cv_master.notify_one();
					}

					if (iNumOfAvailThreads) {
						task.first(iNumOfAvailThreads);
						{
							std::unique_lock<std::mutex> lock(this->taskMutex);
							this->availableCores.push(iNumOfAvailThreads);
						}
					}
					else {
						task.first(task.second);
						{
							std::unique_lock<std::mutex> lock(this->taskMutex);
							this->availableCores.push(task.second);
						}
					}
				}
			}
		)
		);
	}
};

#include <list>

class WorkerQueueWithIDs {
	vector<unique_ptr<thread> > pool; // this holds the threads
	queue< function<void(const int32_t&)> > tasks; // this holds the tasks
	list<int32_t> IDs; // this holds the available IDs
	condition_variable cv_worker;
	condition_variable cv_master;
	mutex taskMutex, idsMutex;
	bool stop = false, bStarted = false;

public:
	/////////////////////////////////////////////////////////
	inline WorkerQueueWithIDs(const int32_t& iNumOfThreads) {
		for (int32_t i = 0; i < iNumOfThreads; ++i) {
			startWorker(i);
			IDs.push_back(i);
		}
	}

	////////////////////////////////////////////////////////
	inline ~WorkerQueueWithIDs() {
		std::unique_lock<std::mutex> lock(taskMutex);
		stop = true;
		cv_worker.notify_all();
		cv_worker.wait(lock, [this] {
			for (size_t i = 0; i < this->pool.size(); ++i) {
				if (this->pool[i] != nullptr) {
					return false;
				}
			}
			return true;
			});
	}

	////////////////////////////////////////////////////////
	inline void pushTask(const function<void(const int32_t&)>& task) {
		tasks.push([task](const int32_t& id) { task(id); });
	}

	////////////////////////////////////////////////////////
	inline void start() {
		if (!tasks.empty()) {
			bStarted = true;
			cv_worker.notify_all();
		}
	}

	////////////////////////////////////////////////////////
	inline void waitUntilFinished() {
		if (!this->tasks.empty()) {
			unique_lock<std::mutex> lock(this->taskMutex);
			this->cv_master.wait(lock,
				[this] { return this->tasks.empty() && (this->IDs.size() == pool.size()); });
		}
		bStarted = false;
	}

private:
	////////////////////////////////////////////////////////
	inline void startWorker(const size_t& index) {
		pool.emplace_back(new thread(
			[this, index]() {
				while (true) {

					function<void(const int32_t&)> task;
					int32_t id = 0;
					{
						std::unique_lock<std::mutex> lock(this->taskMutex);
						
						this->cv_worker.wait(lock, [this]() {return (!(this->tasks.empty()) && bStarted) || this->stop; });


						if (this->stop) {
							this->pool[index]->detach();
							this->pool[index].reset();
							this->cv_worker.notify_one();
							return;
						}

						if (!this->tasks.empty()) {
							task = this->tasks.front();
							this->tasks.pop();

							std::unique_lock<std::mutex> lock2(this->idsMutex);

							id = this->IDs.front();
							this->IDs.pop_front();
						}
						else {
							continue;
						}
					}

					task(id);
					{
						std::unique_lock<std::mutex> lock2(this->idsMutex);
						this->IDs.push_back(id);
					}

					if (this->tasks.empty() && (this->IDs.size() == pool.size())) {
						//std::unique_lock<std::mutex> lock(this->taskMutex);
						cv_master.notify_one();
					}
				}
			}
		)
		);
	}
};