all: ckbDynamicNodeSet ckbDynamicNodeSetSlotDebug ckbDynamicNodeSetSlot test_bucket_sampling

ckbDynamicNodeSet : ckbDynamicNodeSet.cpp PowerlawDegreeSequence.cpp PowerlawDegreeSequence.h
		    g++ -std=c++11 ckbDynamicNodeSet.cpp PowerlawDegreeSequence.cpp -o ckbDynamicNodeSet
ckbDynamicNodeSetSlotDebug : ckbDynamicNodeSetSlot.cpp PowerlawDegreeSequence.cpp PowerlawDegreeSequence.h BucketSampling.h BucketSampling.cpp
		    g++ -fsanitize=address -fno-omit-frame-pointer -g -std=c++11 ckbDynamicNodeSetSlot.cpp PowerlawDegreeSequence.cpp BucketSampling.cpp -o ckbDynamicNodeSetSlotDebug

ckbDynamicNodeSetSlot : ckbDynamicNodeSetSlot.cpp PowerlawDegreeSequence.cpp PowerlawDegreeSequence.h BucketSampling.h BucketSampling.cpp
		    g++ -DNDEBUG -O3 -std=c++11 ckbDynamicNodeSetSlot.cpp PowerlawDegreeSequence.cpp BucketSampling.cpp -o ckbDynamicNodeSetSlot

test_bucket_sampling : test_bucket_sampling.cpp PowerlawDegreeSequence.cpp PowerlawDegreeSequence.h BucketSampling.h BucketSampling.cpp
		    g++ -std=c++11 test_bucket_sampling.cpp BucketSampling.cpp PowerlawDegreeSequence.cpp -g -o test_bucket_sampling

clean :
	rm -f ckbDynamicNodeSet
	rm -f ckbDynamicNodeSetSlotDebug
	rm -f ckbDynamicNodeSetSlot
	rm -f test_bucket_sampling
