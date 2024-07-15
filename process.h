#include <list>

class Process {
public:
    char name;
    std::list<int> cpu_burst_time;
    std::list<int> io_burst_time;
    int next_event_time;
    int cpu_time_left;
    int tau;
    int ready_time;
};

bool operator<(const Process& p1, const Process& p2) {
    return p1.next_event_time > p2.next_event_time ||
        (p1.next_event_time == p2.next_event_time && p1.name > p2.name);
}
