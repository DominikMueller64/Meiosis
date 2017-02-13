#include <random>
#include <boost/uuid/uuid.hpp>            // uuid class
#include <boost/uuid/uuid_generators.hpp> // generators

#ifndef GLOBALS_H
#define GLOBALS_H

/* extern const unsigned random_seed; */
extern std::mt19937 engine;

extern boost::uuids::random_generator uuid_gen;

#endif
