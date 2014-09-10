#include "CS207/Util.hpp"
#include "math.h"
#include <list>

/** Return true iff @a n is prime.
 * @pre @a n >= 0
 */
bool is_prime(int n)
{
  assert(n >= 0);

  static std::list<int> list_of_primes;

  double square_root = sqrt(n);
  


  for (std::list<int>::iterator num = list_of_primes.begin(); *num <= square_root &&  num!=list_of_primes.end(); ++num){
    if (n % *num == 0)
      return false;
  }


  list_of_primes.push_back(n);
  return true;
}

int main()
{
  while (!std::cin.eof()) {
    // How many primes to test? And should we print them?
    std::cerr << "Input Number: ";
    int n = 0;
    CS207::getline_parsed(std::cin, n);
    if (n <= 0)
      break;

    std::cerr << "Print Primes (y/n): ";
    char confirm = 'n';
    CS207::getline_parsed(std::cin, confirm);
    bool print_primes = (confirm == 'y' || confirm == 'Y');

    CS207::Clock timer;


    // Loop and count primes from 2 up to n
    int num_primes = 0;
    if (print_primes) {
      for (int i = 2; i <= n; ++i) {
        if (is_prime(i)) {
          // add i to stack or some data structure
          ++num_primes;
          std::cout << i << std::endl;
       }
      }
    }

    double elapsed_time = timer.seconds();

    std::cout << "There are " << num_primes
              << " primes less than or equal to " << n << ".\n"
              << "Found in " << (1000 * elapsed_time) << " milliseconds.\n\n";
  }

  return 0;
}
