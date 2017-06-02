#ifndef SNPLIB_BARRIER_H
#define SNPLIB_BARRIER_H
#include <atomic>
#include <condition_variable>
#include <mutex>

class Barrier {
private:
  uint32_t num_threads_;
  uint32_t spaces_;
  uint32_t generation_;
  bool is_main_done_;
  std::unique_lock<std::mutex> lk_main_;
  std::mutex mx_main_, mx_threads_;
  std::condition_variable cv_main_, cv_threads_;
  static Barrier *instance_;

  Barrier() = default;

public:
  ~Barrier() = default;
  static Barrier *get_instance() noexcept;
  void SetThreads(const uint32_t num_threads) noexcept;
  void LockMain();
  void ReleaseMain();
  void WaitMain();
  void WaitThreads();
};
#endif
