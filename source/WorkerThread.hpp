/***************************************************************************
*  Part of kASA: https://github.com/SilvioWeging/kASA
*  Inspired by https://github.com/log4cplus/ThreadPool/blob/master/ThreadPool.h
*
*  Copyright (C) 2019 Silvio Weging <silvio.weging@gmail.com>
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

	vector<thread> T;
	int32_t ID = 0;
	vector<function<void(const int32_t&)>> tasks;
	condition_variable cv_worker;
	condition_variable cv_master;
	mutex taskMutex;
	bool stop = false;

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
		cv_worker.notify_one();
	}

	////////////////////////////////////////////////////////
	inline void setNumberOfTasks(const size_t& number) {
		tasks.reserve(number);
	}

	////////////////////////////////////////////////////////
	inline void waitUntilFinished() {
		unique_lock<std::mutex> lock(this->taskMutex);
		this->cv_master.wait(lock,
			[this] { return this->tasks.empty(); });
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
					cv_master.notify_one();
				}
			}
		);
	}

};