#include <gtest/gtest.h>
#include <iostream>
#include <main.h>

#define N 6


namespace {
	class Test_s_cutter : public ::testing::Test {
	 protected:

	  Test_s_cutter() { 
        vector<int> f (N, 2);
        vector<vector<int> > adjacency_list1 (N);
        for (int i = 0 ; i < N; i ++) {
        	for (int j = i + 1; j < N; j ++) {
        		if (rand() % 10 < 4) {
        		  adjacency_list1[i].push_back(j);
        		  adjacency_list1[j].push_back(i);
        		 }
        	}
        }
	    scutter = new S_cutter(adjacency_list, f);
	  }

	  virtual ~Test_s_cutter() {
	    // You can do clean-up work that doesn't throw exceptions here.
	  }

	  // If the constructor and destructor are not enough for setting up
	  // and cleaning up each test, you can define the following methods:

	  virtual void SetUp() {
	    // Code here will be called immediately after the constructor (right
	    // before each test).
	  }

	  virtual void TearDown() {
	    // Code here will be called immediately after each test (right
	    // before the destructor).
	  }

	  // Objects declared here can be used by all tests in the test case for Foo.
	  S_cutter* scutter;
	};

	// Tests that the Foo::Bar() method does Abc.
	TEST_F(FooTest, MethodBarDoesAbc) {
	  const string input_filepath = "this/package/testdata/myinputfile.dat";
	  const string output_filepath = "this/package/testdata/myoutputfile.dat";
	  Foo f;
	  EXPECT_EQ(0, f.Bar(input_filepath, output_filepath));
	}

	// Tests that Foo does Xyz.
	TEST_F(FooTest, DoesXyz) {
	  // Exercises the Xyz feature of Foo.
} // namespace


int main (int argc, char* argv[]) {
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}