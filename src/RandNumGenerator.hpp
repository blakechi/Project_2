#include <random>

class RandomNumberBetween
{
public:
    RandomNumberBetween(int low, int high)
    : _random_engine{ std::random_device{}() }, _distribution{ low, high }
    {
    }

    int operator()()
    {
        return _distribution(_random_engine);
    }

private:
    std::mt19937 _random_engine;
    std::uniform_int_distribution<int> _distribution;
};


// Examples:
// std::generate(begin(range), end(range), generator);
// std::generate_n(std::back_inserter(numbers), numBlank, RandomNumberBetween(0, dataCount));