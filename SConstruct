env = Environment()

env.Append(CCFLAGS=['-std=c++11'])
env.Append(CCFLAGS=['-O3'])
env.Append(CCFLAGS=['-g'])
env.Append(CCFLAGS=['-Wall'])
#env.Append(CCFLAGS=['-Werror'])
env.Append(CCFLAGS=['-pg'])
env.Append(LINKFLAGS=['-pg'])

main = env.Program("Random",["Random.cpp", "DP.cc"])
env.Program("Random50_45mix", ["Random50_45mix.cpp", "DP.cc"])
