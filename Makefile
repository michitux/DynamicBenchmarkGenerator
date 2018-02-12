all: ckbDynamicNodeSet test_bucket_sampling

ckbDynamicNodeSet : ckbDynamicNodeSet.cpp PowerlawDegreeSequence.cpp PowerlawDegreeSequence.h
		    g++ -std=c++11 ckbDynamicNodeSet.cpp PowerlawDegreeSequence.cpp -o ckbDynamicNodeSet
test_bucket_sampling : test_bucket_sampling.cpp PowerlawDegreeSequence.cpp PowerlawDegreeSequence.h BucketSampling.h BucketSampling.cpp
		    g++ -std=c++11 test_bucket_sampling.cpp BucketSampling.cpp PowerlawDegreeSequence.cpp -g -o test_bucket_sampling

clean :
	rm -f ckbDynamicNodeSet
	rm -f test_bucket_sampling
