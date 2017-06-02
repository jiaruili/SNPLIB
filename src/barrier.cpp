#include "barrier.h"

Barrier *Barrier::instance_ = nullptr;

Barrier *Barrier::get_instance() noexcept {
  if (instance_ == nullptr) {
    instance_ = new Barrier;
  }
  return instance_;
}

void Barrier::SetThreads(const uint32_t num_threads) noexcept {
  num_threads_ = num_threads + 1;
  spaces_ = num_threads + 1;
  generation_ = 0;
  is_main_done_ = false;
  lk_main_ = std::unique_lock<std::mutex>(mx_main_, std::defer_lock);
}

void Barrier::LockMain() { lk_main_.lock(); }

void Barrier::ReleaseMain() {
  is_main_done_ = true;
  lk_main_.unlock();
  cv_main_.notify_all();
}

void Barrier::WaitMain() {
  std::unique_lock<std::mutex> lk(mx_main_);
  while (!is_main_done_) {
    cv_main_.wait(lk);
  }
  lk.unlock();
}

void Barrier::WaitThreads() {
  std::unique_lock<std::mutex> lk(mx_threads_);
  uint32_t local_gen = generation_;
  if (!--spaces_) {
    spaces_ = num_threads_;
    ++generation_;
    is_main_done_ = false;
    lk.unlock();
    cv_threads_.notify_all();
  } else {
    while (local_gen == generation_) {
      cv_threads_.wait(lk);
    }
    lk.unlock();
  }
}
