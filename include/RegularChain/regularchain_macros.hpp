

#include "../Utils/Parallel/AsyncGeneratorThread.hpp"
#include "../Utils/Parallel/AsyncGeneratorPool.hpp"


#if defined(SERIAL) && SERIAL
#define TRIANGULARIZE_SERIAL 1
#else 
#define TRIANGULARIZE_SERIAL 0
#endif


#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS

	//consumer
	#define RC_TRI_CONSUMER_INIT(retType, retVar, methodName, arg1, ...) AsyncGeneratorThread<retType> retVar(methodName, (arg1),  __VA_ARGS__)

	//consumer
	#define RC_TRI_CONSUMER_DECLARE_LOOPELEM(elemType, elemVar) elemType elemVar; 

	//consumer
	#define RC_TRI_CONSUMER_LOOP(loopVec, loopElem, i) {size_t i = 0;} while(loopVec.getNextObject((loopElem)))

	//consumer
	#define RC_TRI_CONSUMER_GET_LOOPELEM(currElem, loopVec, i)

	//consumer
	#define RC_TRI_COMSUMER_IS_LAST_ITER(loopVec, i) ((loopVec).isFinished())


	//producer
	#define RC_TRI_PRODUCER_DECLARE_RESULT(retType, result) 

	//producer
	#define RC_TRI_PRODUCER_ACCUMULATE(accumulator, objToAcc) ((accumulator).generateObject((objToAcc)))

	//producer
	#define RC_TRI_PRODUCER_COMPLETE(result) result.setComplete(); return

#else

	//consumer
	#define RC_TRI_CONSUMER_INIT(retType, retVar, methodName, arg1, ...) std::vector<retType> retVar = ((*arg1).*(methodName))( __VA_ARGS__ )

	//consumer
	#define RC_TRI_CONSUMER_DECLARE_LOOPELEM(elemType, elemVar) elemType elemVar; 

	//consumer
	#define RC_TRI_CONSUMER_LOOP(loopVec, loopElem, i) for(size_t (i) = 0; (i) < (loopVec).size(); ++(i))

	//consumer	
	#define RC_TRI_CONSUMER_GET_LOOPELEM(loopVec, currElem, i) ((currElem) = (loopVec)[(i)]) 

	//consumer
	#define RC_TRI_COMSUMER_IS_LAST_ITER(loopVec, i) ((i) + 1 == (loopVec).size())


	//producer
	#define RC_TRI_PRODUCER_DECLARE_RESULT(retType, result) std::vector<retType> result; 

	//producer
	#define RC_TRI_PRODUCER_ACCUMULATE(accumulator, objToAcc) ((accumulator).push_back((objToAcc)))

	//producer
	#define RC_TRI_PRODUCER_COMPLETE(result) return result

#endif



#if defined(RC_WITH_GENERATORS) && RC_WITH_GENERATORS

	//consumer
	#define RC_GEN_CONSUMER_INIT(retType, retVar, methodName, arg1, ...) AsyncGeneratorPool<retType> retVar(methodName, (arg1),  __VA_ARGS__)

	//consumer
	#define RC_GEN_CONSUMER_DECLARE_LOOPELEM(elemType, elemVar) elemType elemVar; 

	//consumer
	#define RC_GEN_CONSUMER_LOOP(loopVec, loopElem, i) {size_t i = 0;} while((loopVec).getNextObject((loopElem)))

	//consumer
	#define RC_GEN_CONSUMER_GET_LOOPELEM(currElem, loopVec, i)

	//consumer
	#define RC_GEN_COMSUMER_IS_LAST_ITER(loopVec, i) ((loopVec).isFinished())


	//producer
	#define RC_GEN_PRODUCER_DECLARE_RESULT(retType, result) 

	//producer
	#define RC_GEN_PRODUCER_ACCUMULATE(accumulator, objToAcc) ((accumulator).generateObject((objToAcc)))

	//producer
	#define RC_GEN_PRODUCER_COMPLETE(result) result.setComplete(); return


#else

	//consumer
	#define RC_GEN_CONSUMER_INIT(retType, retVar, methodName, arg1, ...) std::vector<retType> retVar = ((*arg1).*(methodName))( __VA_ARGS__ )
	#define RC_GEN_FORCE_CONSUMER_INIT(retType, retVar, methodName, arg1, ...) std::vector<retType> retVar = ((*arg1).*(methodName))( __VA_ARGS__ )

	//consumer
	#define RC_GEN_CONSUMER_DECLARE_LOOPELEM(elemType, elemVar) elemType elemVar; 

	//consumer
	#define RC_GEN_CONSUMER_LOOP(loopVec, loopElem, i) for(size_t (i) = 0; (i) < (loopVec).size(); ++(i))

	//consumer	
	#define RC_GEN_CONSUMER_GET_LOOPELEM(loopVec, currElem, i) ((currElem) = (loopVec)[(i)]) 

	//consumer
	#define RC_GEN_COMSUMER_IS_LAST_ITER(loopVec, i) ((i) + 1 == (loopVec).size())

	//producer
	#define RC_GEN_PRODUCER_DECLARE_RESULT(retType, result) std::vector<retType> result; 

	//producer
	#define RC_GEN_PRODUCER_ACCUMULATE(accumulator, objToAcc) ((accumulator).push_back((objToAcc)))

	//producer
	#define RC_GEN_PRODUCER_COMPLETE(result) return result


#endif
