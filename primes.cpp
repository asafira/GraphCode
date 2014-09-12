#include "CS207/Util.hpp"
#include "math.h"
#include <list>

/** @brief Return true iff @a n is prime.
 *  @pre @a n >= 0
 *  @post true if there does not exist 2 <= x <= n s.t. x | n for int x, 
 *  false otherwise.
 */
bool is_prime(int n)
{
  assert(n >= 0); // Check for a non-negative integer

  // Keep  persistent prime list
  static std::list<int> list_of_primes;

  // Calculute square_root for better look bounds
  int square_root = (int) floor(sqrt(n)); 

  // Look through list of primes until it ends or we reach the sq. root
  for (std::list<int>::iterator num = list_of_primes.begin(); *num <= square_root &&  num != list_of_primes.end(); ++num){
    
    // Return false if we find a factor
    if (n % *num == 0)
      return false;
  }


  // Add to prime list if we cannot find a factor, return true
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
