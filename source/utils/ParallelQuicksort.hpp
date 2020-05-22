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


///////////////////////////////////////////////////////
// The following code is a slightly modified version of the working draft from Microsoft
// for the parallel stl implementation
// see https://archive.codeplex.com/?p=parallelstl for further details
// It was published under the Apache License 2.0
//
//	Copyright 2020 Silvio Weging
//
//	Licensed under the Apache License, Version 2.0 (the "License");
//	you may not use this file except in compliance with the License.
//	You may obtain a copy of the License at
//
//	https://www.apache.org/licenses/LICENSE-2.0
//
//	Unless required by applicable law or agreed to in writing, software
//	distributed under the License is distributed on an "AS IS" BASIS,
//	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//	See the License for the specific language governing permissions and
//	limitations under the License.


#pragma once

#include "../MetaHeader.h"

namespace Utilities {

	//////////////////////////////////////////////////////////////////////////////////
	// This number is used to control dynamic task splitting
	// The ideal chunk (task) division is that the number of cores is equal to the number of tasks, but it will 
	// perform very poorly when tasks are not balanced. The simple solution is to allocate more tasks than number 
	// of cores. _SortMaxTasksPerCore provides a maximum number of tasks that will be allocated per core.
	// If this number is too small, the load balancing problem will affect efficiency very seriously, especially
	// when the compare operation is expensive.
	//
	// Note that this number is a maximum number -- the dynamic partition system will reduce the number of partitions
	// per core based on the dynamic load. If all cores are very busy, the number of partitions will shrink to 
	// reduce the scheduler overhead.
	//
	// Initially, the total tasks(chunks) number of partitions "_Div_num" will be: core number * _SortMaxTasksPerCore.
	// The _Div_num will be divided by 2 after each task splitting. There are two special numbers for _Div_num:
	//     1. When _Div_num reaches the point that _Div_num < _SortMaxTasksPerCore, it means we have split more tasks than cores.
	//     2. When _Div_num reaches the point that _Div_num <= 1, it means stop splitting more tasks and begin sorting serially.
	const int _SortMaxTasksPerCore = 8; // original: 1024

	// This is a number mainly is used to control the sampling and dynamic task splitting strategies.
	// If the user configurable minimal divisible chunk size (default is 2048) is smaller than FINE_GRAIN_CHUNK_SIZE,
	// the random sampling algorithm for quicksort will enter fine-grained mode, and take a strategy that reduces the sampling 
	// overhead. Also, the dynamic task splitting will enter fine-grained mode, which will split as many tasks as possible.
	const int _SortChunkSize = 512;

	// This is the maximum depth that the quicksort will be called recursively. If we allow too far, a stack overflow may occur.
	const int _SortMaxRecursionDepth = 8; // original: 64

	// Comment from Silvio Weging: 
	// increasing the _SortMaxTasksPerCore and _SortMaxRecursionDepth results in a lot of virtual memory allocation (without real memory allocations).
	// This effect is visible in linux since the stack allocated for a thread is 8MB and a lot of threads are getting created without any real usage of that memory.
	// Since kASA assumes an x64 system, this -should- cause no issues.


	///////////////////////////////////////////////////////
	template<typename _Random_iterator, typename _Function>
	inline size_t _Median_of_three(const _Random_iterator &_Begin, const size_t& _A, const size_t& _B, const size_t& _C, _Function &_Func, bool &_Potentially_equal)
	{
		_Potentially_equal = false;
		if (_Func(_Begin[_A], _Begin[_B]))
		{
			if (_Func(_Begin[_A], _Begin[_C]))
			{
				return _Func(_Begin[_B], _Begin[_C]) ? _B : _C;
			}
			else
			{
				return _A;
			}
		}
		else
		{
			if (_Func(_Begin[_B], _Begin[_C]))
			{
				return _Func(_Begin[_A], _Begin[_C]) ? _A : _C;
			}
			else
			{
				_Potentially_equal = true;
				return _B;
			}
		}
	}

	///////////////////////////////////////////////////////

	template<typename _Random_iterator, typename _Function>
	inline size_t _Median_of_nine(const _Random_iterator &_Begin, const size_t& _Size, _Function &_Func, bool &_Potentially_equal)
	{
		size_t _Offset = _Size / 8;
		size_t _A = _Median_of_three(_Begin, 0, _Offset, _Offset * 2, _Func, _Potentially_equal),
			_B = _Median_of_three(_Begin, _Offset * 3, _Offset * 4, _Offset * 5, _Func, _Potentially_equal),
			_C = _Median_of_three(_Begin, _Offset * 6, _Offset * 7, _Size - 1, _Func, _Potentially_equal);
		_B = _Median_of_three(_Begin, _A, _B, _C, _Func, _Potentially_equal);

		if (_Potentially_equal)
		{
			_Potentially_equal = !_Func(_Begin[_C], _Begin[_A]);
		}

		return _B;
	}

	///////////////////////////////////////////////////////

	// _Potentially_equal means that potentially all the values in the buffer are equal to the pivot value
	template<typename _Random_iterator, typename _Function>
	inline size_t _Select_median_pivot(const _Random_iterator &_Begin, const size_t& _Size, _Function &_Func, const size_t& _Chunk_size, bool &_Potentially_equal)
	{
		// Base on different chunk size, apply different sampling optimization
		if (_Chunk_size < _SortChunkSize && _Size <= std::max<size_t>(_Chunk_size * 4, static_cast<size_t>(15)))
		{
			bool _Never_care_equal;
			return _Median_of_three(_Begin, 0, _Size / 2, _Size - 1, _Func, _Never_care_equal);
		}
		else
		{
			return _Median_of_nine(_Begin, _Size, _Func, _Potentially_equal);
		}
	}


	///////////////////////////////////////////////////////
	template<typename _Random_iterator, typename _Function>
	void _Parallel_quicksort_impl(const _Random_iterator &_Begin, const size_t& _Size, _Function &_Func, const size_t& _Div_num, const size_t& _Chunk_size, const int& _Depth)
	{
		// Modification: Added some paranthesis. 2020 Silvio Weging
		if ((_Depth >= _SortMaxRecursionDepth) || (_Size <= _Chunk_size) || (_Size <= static_cast<size_t>(3)) || ((_Chunk_size >= _SortChunkSize) && (_Div_num <= 1)))
		{
			return std::sort(_Begin, _Begin + _Size, _Func);
		}

		// Determine whether we need to do a three-way quick sort
		// We benefit from three-way merge if there are a lot of elements that are EQUAL to the median value,
		// _Select_median_pivot function will test redundant density by sampling
		bool _Is_three_way_split = false;
		size_t _Mid_index = _Select_median_pivot(_Begin, _Size, _Func, _Chunk_size, _Is_three_way_split);

		// Move the median value to the _Begin position.
		if (_Mid_index)
		{
			std::swap(*_Begin, _Begin[_Mid_index]);
		}
		size_t _I = 1, _J = _Size - 1;

		// Three-way or two-way partition
		// _Div_num < _SortMaxTasksPerCore is checked to make sure it will never do three-way split before splitting enough tasks
		if (_Is_three_way_split && _Div_num < _SortMaxTasksPerCore)
		{
			while (_Func(*_Begin, _Begin[_J]))
			{
				--_J;
			}

			while (_Func(_Begin[_I], *_Begin))
			{
				++_I;
			}

			// Starting from this point, left side of _I will less than median value, right side of _J will be greater than median value, 
			// and the middle part will be equal to median. _K is used to scan between _I and _J
			size_t _K = _J;
			while (_I <= _K)
			{
				if (_Func(_Begin[_K], *_Begin))
				{
					std::swap(_Begin[_I++], _Begin[_K]);
				}
				else
				{
					--_K;
				}

				while (_Func(*_Begin, _Begin[_K]))
				{
					std::swap(_Begin[_K--], _Begin[_J--]);
				}
			}

			++_J;
		}
		else
		{
			while (_I <= _J)
			{
				// Will stop before _Begin
				while (_Func(*_Begin, _Begin[_J]))
				{
					--_J;
				}

				// There must be another element equal or greater than *_Begin
				while (_Func(_Begin[_I], *_Begin))
				{
					++_I;
				}

				if (_I < _J)
				{
					std::swap(_Begin[_I++], _Begin[_J--]);
				}
				else
				{
					break;
				}
			}

			_I = ++_J;
		}

		std::swap(*_Begin, _Begin[--_I]);

		//////////////////////////////////////////////////
		// Further modifications start here... (2020, Silvio Weging)

		size_t _Next_div = _Div_num / 2;
		auto handle = ([&]
		{
			_Parallel_quicksort_impl(_Begin + _J, _Size - _J, _Func, _Next_div, _Chunk_size, _Depth + 1);
		});

		thread t(handle);

		_Parallel_quicksort_impl(_Begin, _I, _Func, _Next_div, _Chunk_size, _Depth + 1);


		// If at this point, the work hasn't been scheduled, then slow down creating new tasks
		if (_Div_num < _SortMaxTasksPerCore)
		{
			_Next_div /= 2;
		}

		t.join();
	}

	template<typename _FwdIt, typename _Pr>
	inline void parallelQuicksort(_FwdIt _First, _FwdIt _Last, _Pr _Pred, const int32_t& iNumOfThreads)
	{
		// Check for cancellation before the algorithm starts.
		size_t _Size = _Last - _First;
		size_t _Core_num = static_cast<size_t>(iNumOfThreads);
		const size_t _ChunkSize = 2048; // Default chunk size

		if (_Size <= _ChunkSize || _Core_num < 2)
		{
			return std::sort(_First, _Last, _Pred);
		}

		_Parallel_quicksort_impl(_First, _Size, _Pred, _Core_num * _SortMaxTasksPerCore, _ChunkSize, 0);
	}
	//////////////////////////////////////////////////
	// ... and end here.
}