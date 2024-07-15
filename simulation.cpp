#include <string>
#include <queue>
#include <list>
#include <iostream>
#include <climits>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include "process.h"



struct tau_time {
    bool operator() (const Process& p1, const Process& p2) {
        return p1.tau > p2.tau || (p1.tau == p2.tau && p1.name > p2.name);
    }
};

struct adjusted_tau_time {
    bool operator() (const Process& p1, const Process& p2) {
        int p1_time_used = p1.cpu_burst_time.front() - p1.cpu_time_left;
        int p2_time_used = p2.cpu_burst_time.front() - p2.cpu_time_left;
        return ((p1.tau - p1_time_used) > (p2.tau - p2_time_used) ||
        ((p1.tau - p1_time_used) == (p2.tau - p2_time_used) && p1.name > p2.name));
    }
};

void print_list_queue(const std::list<Process>& ready_queue) {
    std::cout << " [Q:";
    if (ready_queue.empty()) {
        std::cout << " empty]" << std::endl;
    }
    else {
        for (std::list<Process>::const_iterator itr = ready_queue.begin(); itr != ready_queue.end(); itr++) {
            std::cout << " " << itr->name;
        }
        std::cout << "]" << std::endl;
    }
}

void print_priority_queue(const std::priority_queue<Process,std::vector<Process>,tau_time>& ready_queue) {
    std::cout << " [Q:";
    if (ready_queue.empty()) {
        std::cout << " empty]" << std::endl;
    }
    else {
        std::priority_queue<Process,std::vector<Process>,tau_time> copy_queue(ready_queue);
        while (!copy_queue.empty()) {
            std::cout << " " << copy_queue.top().name;
            copy_queue.pop();
        }
        std::cout << "]" << std::endl;
    }
}

void print_priority_queue(const std::priority_queue<Process,std::vector<Process>,adjusted_tau_time>& ready_queue) {
    std::cout << " [Q:";
    if (ready_queue.empty()) {
        std::cout << " empty]" << std::endl;
    }
    else {
        std::priority_queue<Process,std::vector<Process>,adjusted_tau_time> copy_queue(ready_queue);
        while (!copy_queue.empty()) {
            std::cout << " " << copy_queue.top().name;
            copy_queue.pop();
        }
        std::cout << "]" << std::endl;
    }
}

void print_time_and_process(int time, char name) {
    std::cout << "time " << time << "ms: Process " << name;
}

void print_algorithm_stats(std::ofstream& output, int total_number_bursts, int wait_time, int turnaround_time, int total_bursts_time,
                            int context_switches, int preemptions, int total_time, int FCFS, int SJF, int SRT, int RR)
{
    output << "Algorithm ";
    if (FCFS == 1) output << "FCFS" << std::endl;
    if (SJF == 1) output << "SJF" << std::endl;
    if (SRT == 1) output << "SRT" << std::endl;
    if (RR == 1) output << "RR" << std::endl;
    float avg_wait = wait_time / float(total_number_bursts);
    float avg_turnaround = turnaround_time / float(total_number_bursts);
    float avg_burst = total_bursts_time / float(total_number_bursts);
    float cpu_utilization = (total_bursts_time / float(total_time)) * 100;
    output << std::fixed << std::setprecision(3);
    avg_burst += 0.000499;
    avg_wait += 0.000499;      //for rounding purposes
    avg_turnaround += 0.000499;
    cpu_utilization += 0.000499;
    output << "-- average CPU burst time: " << avg_burst << " ms" << std::endl;
    output << "-- average wait time: " << avg_wait << " ms" << std::endl;
    output << "-- average turnaround time: " << avg_turnaround  << " ms" << std::endl;
    output << "-- total number of context switches: " << context_switches << std::endl;
    output << "-- total number of preemptions: " << preemptions << std::endl;
    output << "-- CPU utilization: " << cpu_utilization << "%" << std::endl;

}

int recalculate_tau(int old_tau, int runtime, double alpha) {
    return ceil(alpha*runtime + (1-alpha)*old_tau);
}

double next_exp(double l, double b)
{
	double r = drand48();
	double x = -log( r ) / l; //x = -ln(r)/lambda

	if(x <= b)
		return x;

	return next_exp(l, b);
}

void fcfs(const std::vector<Process>& process_list, int context_time, std::ofstream& output) {
    Process* running = nullptr;
    std::list<Process> ready_queue;
    std::priority_queue<Process> io_wait;
    std::priority_queue<Process> unstarted;
    int cpu_ready_time = context_time/2;
    int next_time = 0;

    std::vector<int> start_times;
    int total_bursts_time = 0;
    int total_number_bursts = 0;
    int wait_time = 0;
    int turnaround_time = 0;
    int context_switches = 0;
    int preemptions = 0;
    int total_time = 0;
    //int current_time_stamp = 0;
    for (unsigned int i = 0; i < process_list.size(); i++) {
        unstarted.push(process_list[i]);
        total_number_bursts += process_list[i].cpu_burst_time.size();
    }
    std::cout << "time 0ms: Simulator started for FCFS [Q: empty]" << std::endl;
    while (running != nullptr || ready_queue.size() > 0 || io_wait.size() > 0 || unstarted.size() > 0) {
        next_time = INT_MAX;
        std::string next_event = "";
        if (running != nullptr) {
            next_time = running->next_event_time;
            next_event = "cpu burst complete";
        }
        else if (!ready_queue.empty()) {
            if (cpu_ready_time < next_time) {
                next_time = cpu_ready_time;
                next_event = "entering cpu burst";
            }
        }
        if (!io_wait.empty()) {
            Process p = io_wait.top();
            if (p.next_event_time < next_time) {
                next_time = p.next_event_time;
                next_event = "io burst complete";
            }
        }
        if (!unstarted.empty()) {
            Process p = unstarted.top();
            if (p.next_event_time < next_time) {
                next_time = p.next_event_time;
                next_event = "process enters the scheduler";
            }
        }

        if (next_event == "cpu burst complete") {
            //std::cout << "next time: " << next_time << std::endl;
            //std::cout << "arrival time: " << arrival_time << std::endl;
            //total_turnaround_time += next_time;
            //turnaround_time += (next_time - arrival_time);
            //wait_time += (next_time - arrival_time) - running->cpu_burst_time.front();
            //std::cout << "turnaround time: " << (next_time - arrival_time) << std::endl;
            //std::cout << "wait time: " << (next_time - arrival_time) - running->cpu_burst_time.front() << std::endl;
            (running->cpu_burst_time).pop_front();

            if ((running->io_burst_time).size() > 0) {
                if (next_time < 1000) {
                    print_time_and_process(next_time,running->name);
                    std::cout << " completed a CPU burst; " << (running->cpu_burst_time).size() << " burst";
                    if ((running->cpu_burst_time).size() != 1) std::cout << "s";
                    std::cout << " to go";
                    print_list_queue(ready_queue);
                    print_time_and_process(next_time,running->name);
                    std::cout << " switching out of CPU; will block on I/O until time " << next_time + context_time/2 + (running->io_burst_time).front() << "ms";
                    print_list_queue(ready_queue);
                }

                turnaround_time += next_time - running->ready_time + context_time/2;
                running->next_event_time = next_time + context_time/2 + (running->io_burst_time).front();
                io_wait.push(*running);
            }
            else {
                turnaround_time += next_time - running->ready_time + context_time/2;

                print_time_and_process(next_time,running->name);
                std::cout << " terminated";
                print_list_queue(ready_queue);
            }

            delete(running);
            running = nullptr;
            cpu_ready_time = next_time + context_time;
        }

        else if (next_event == "entering cpu burst") {
            running = new Process();
            *running = ready_queue.front();
            running->next_event_time = next_time + (running->cpu_burst_time).front();
            ready_queue.pop_front();

            wait_time += next_time - context_time / 2;
            total_bursts_time += running->cpu_burst_time.front();
            context_switches++;
            if (next_time < 1000) {
                print_time_and_process(next_time,running->name);
                std::cout << " started using the CPU for " << (running->cpu_burst_time).front() << "ms burst";
                print_list_queue(ready_queue);
            }
        }

        else if (next_event == "io burst complete") {
            Process p = io_wait.top();
            io_wait.pop();
            p.io_burst_time.pop_front();
            if (ready_queue.empty() && cpu_ready_time < next_time + context_time/2) {
                cpu_ready_time = next_time + context_time/2;
            }
            p.ready_time = next_time;
            wait_time -= next_time;
            ready_queue.push_back(p);

            if (next_time < 1000) {
                print_time_and_process(next_time,p.name);
                std::cout << " completed I/O; added to ready queue";
                print_list_queue(ready_queue);
            }
            start_times.push_back(next_time);
        }

        else if (next_event == "process enters the scheduler") {
            Process p = unstarted.top();
            unstarted.pop();
            if (ready_queue.empty() && cpu_ready_time < next_time + context_time/2) {
                cpu_ready_time = next_time + context_time/2;
            }
            wait_time -= next_time;
            p.ready_time = next_time;
            ready_queue.push_back(p);
            if (next_time < 1000) {
                print_time_and_process(next_time,p.name);
                std::cout << " arrived; added to ready queue";
                print_list_queue(ready_queue);
            }
        }
        else {
            std::cout << "Unknown event occurred." << std::endl;
        }
    }
    total_time = next_time + context_time/2;
    std::cout << "time " << next_time + context_time/2 << "ms: Simulator ended for FCFS [Q: empty]" << std::endl;
    print_algorithm_stats(output,total_number_bursts,wait_time,turnaround_time,total_bursts_time,
                            context_switches,preemptions,total_time,1,0,0,0);
}

void sjf(const std::vector<Process>& process_list, int context_time, double alpha, std::ofstream& output) {
    Process* running = nullptr;
    std::priority_queue<Process,std::vector<Process>,tau_time> ready_queue;
    std::priority_queue<Process> io_wait;
    std::priority_queue<Process> unstarted;
    int cpu_ready_time = context_time/2;
    int next_time = 0;
    int choose_next_process = 0;

    int total_bursts_time = 0;
    int total_number_bursts = 0;
    int wait_time = 0;
    //int arrival_time;
    //int total_turnaround_time = 0;
    int turnaround_time = 0;
    int context_switches = 0;
    int preemptions = 0;
    int total_time = 0;
    for (unsigned int i = 0; i < process_list.size(); i++) {
        unstarted.push(process_list[i]);
        total_number_bursts += process_list[i].cpu_burst_time.size();
    }
    std::cout << "time 0ms: Simulator started for SJF [Q: empty]" << std::endl;
    while (running != nullptr || ready_queue.size() > 0 || io_wait.size() > 0 || unstarted.size() > 0) {
        next_time = INT_MAX;
        std::string next_event = "";
        if (running != nullptr && choose_next_process == INT_MAX) {
            next_time = running->next_event_time;
            next_event = "cpu burst complete";
        }
        else if (!ready_queue.empty() || (choose_next_process != INT_MAX && running != nullptr)) {
            if (cpu_ready_time < next_time) {
                next_time = cpu_ready_time;
                next_event = "entering cpu burst";
            }
        }
        if (!io_wait.empty()) {
            Process p = io_wait.top();
            if (p.next_event_time < next_time) {
                next_time = p.next_event_time;
                next_event = "io burst complete";
            }
        }
        if (!unstarted.empty()) {
            Process p = unstarted.top();
            if (p.next_event_time < next_time) {
                next_time = p.next_event_time;
                next_event = "process enters the scheduler";
            }
        }

        if (next_time > choose_next_process && !ready_queue.empty() && running == nullptr) {
            running = new Process();
            *running = ready_queue.top();
            ready_queue.pop();
        }

        if (next_event == "cpu burst complete") {
            int runtime = (running->cpu_burst_time).front();
            //std::cout << "next time: " << next_time << std::endl;
            //std::cout << "arrival time: " << arrival_time << std::endl;
            //total_turnaround_time += next_time;
            //turnaround_time += (next_time - arrival_time);
            //wait_time += (next_time - arrival_time) - running->cpu_burst_time.front();
            //std::cout << "turnaround time: " << (next_time - arrival_time) << std::endl;
            //std::cout << "wait time: " << (next_time - arrival_time) - running->cpu_burst_time.front() << std::endl;

            (running->cpu_burst_time).pop_front();

            if ((running->io_burst_time).size() > 0) {
                if (next_time < 1000) {
                    print_time_and_process(next_time,running->name);
                    std::cout << " (tau " << running->tau << "ms) completed a CPU burst; " << (running->cpu_burst_time).size() << " burst";
                    if ((running->cpu_burst_time).size() != 1) std::cout << "s";
                    std::cout << " to go";
                    print_priority_queue(ready_queue);
                }

                int new_tau = recalculate_tau(running->tau,runtime,alpha);
                if (next_time < 1000) {
                    std::cout << "time " << next_time << "ms: Recalculated tau for process " << running->name << ": old tau " << running->tau << "ms; new tau " << new_tau << "ms";
                    print_priority_queue(ready_queue);
                }
                running->tau = new_tau;

                if (next_time < 1000) {
                    print_time_and_process(next_time,running->name);
                    std::cout << " switching out of CPU; will block on I/O until time " << next_time + context_time/2 + (running->io_burst_time).front() << "ms";
                    print_priority_queue(ready_queue);
                }
                turnaround_time += next_time - running->ready_time + context_time/2;
                running->next_event_time = next_time + context_time/2 + (running->io_burst_time).front();
                io_wait.push(*running);

            }
            else {
                print_time_and_process(next_time,running->name);
                std::cout << " terminated";
                print_priority_queue(ready_queue);

                turnaround_time += next_time - running->ready_time + context_time/2;
            }

            delete(running);
            running = nullptr;
            cpu_ready_time = next_time + context_time;
            choose_next_process = next_time + context_time/2;
        }

        else if (next_event == "entering cpu burst") {
            running->next_event_time = next_time + (running->cpu_burst_time).front();
            choose_next_process = INT_MAX;

            wait_time += next_time - context_time/2;
            total_bursts_time += running->cpu_burst_time.front();
            context_switches++;
            if (next_time < 1000) {
                print_time_and_process(next_time,running->name);
                std::cout << " (tau " << running->tau << "ms) started using the CPU for " << (running->cpu_burst_time).front() << "ms burst";
                print_priority_queue(ready_queue);
            }
        }

        else if (next_event == "io burst complete") {
            Process p = io_wait.top();
            io_wait.pop();
            p.io_burst_time.pop_front();
            if (ready_queue.empty() && running == nullptr && cpu_ready_time < next_time + context_time/2) {
                choose_next_process = next_time;
                cpu_ready_time = next_time + context_time/2;
            }
            p.ready_time = next_time;
            wait_time -= next_time;
            ready_queue.push(p);

            if (next_time < 1000) {
                print_time_and_process(next_time,p.name);
                std::cout << " (tau " << p.tau << "ms) completed I/O; added to ready queue";
                print_priority_queue(ready_queue);
            }
        }

        else if (next_event == "process enters the scheduler") {
            Process p = unstarted.top();
            unstarted.pop();
            if (ready_queue.empty() && running == nullptr && cpu_ready_time < next_time + context_time/2) {
                choose_next_process = next_time;
                cpu_ready_time = next_time + context_time/2;
            }
            p.ready_time = next_time;
            wait_time -= next_time;
            ready_queue.push(p);

            if (next_time < 1000) {
                print_time_and_process(next_time,p.name);
                std::cout << " (tau " << p.tau << "ms) arrived; added to ready queue";
                print_priority_queue(ready_queue);
            }
        }
        else {
            std::cout << "Unknown event occurred." << std::endl;
        }
    }
    total_time = next_time + context_time / 2;
    std::cout << "time " << next_time + context_time/2 << "ms: Simulator ended for SJF [Q: empty]" << std::endl;
    print_algorithm_stats(output,total_number_bursts,wait_time,turnaround_time,total_bursts_time,
                            context_switches,preemptions,total_time,0,1,0,0);
}

void srt(const std::vector<Process>& process_list, int context_time, double alpha, std::ofstream& output) {
    Process* running = nullptr;
    std::priority_queue<Process,std::vector<Process>,adjusted_tau_time> ready_queue;
    std::priority_queue<Process> io_wait;
    std::priority_queue<Process> unstarted;
    int cpu_ready_time = context_time/2;
    int next_time = 0;
    int choose_next_process = 0;

    int total_bursts_time = 0;
    int total_number_bursts = 0;
    int wait_time = 0;
    int turnaround_time = 0;
    int context_switches = 0;
    int preemptions = 0;
    int total_time = 0;
    for (unsigned int i = 0; i < process_list.size(); i++) {
        unstarted.push(process_list[i]);
        total_number_bursts += process_list[i].cpu_burst_time.size();
    }
    std::cout << "time 0ms: Simulator started for SRT [Q: empty]" << std::endl;
    while (running != nullptr || ready_queue.size() > 0 || io_wait.size() > 0 || unstarted.size() > 0) {
        next_time = INT_MAX;
        std::string next_event = "";
        if (running != nullptr && choose_next_process == INT_MAX) {
            next_time = running->next_event_time;
            next_event = "cpu burst complete";
        }
        else if (!ready_queue.empty() || (choose_next_process != INT_MAX && running != nullptr)) {
            if (cpu_ready_time < next_time) {
                next_time = cpu_ready_time;
                next_event = "entering cpu burst";
            }
        }
        if (!io_wait.empty()) {
            Process p = io_wait.top();
            if (p.next_event_time < next_time) {
                next_time = p.next_event_time;
                next_event = "io burst complete";
            }
        }
        if (!unstarted.empty()) {
            Process p = unstarted.top();
            if (p.next_event_time < next_time) {
                next_time = p.next_event_time;
                next_event = "process enters the scheduler";
            }
        }

        if (next_time > choose_next_process && !ready_queue.empty() && running == nullptr) {
            running = new Process();
            *running = ready_queue.top();
            ready_queue.pop();
        }

        if (next_event == "cpu burst complete") {
            turnaround_time += next_time + context_time/2;
            int runtime = (running->cpu_burst_time).front();
            total_bursts_time += runtime;
            (running->cpu_burst_time).pop_front();

            if ((running->io_burst_time).size() > 0) {
                running->cpu_time_left = (running->cpu_burst_time).front();
                if (next_time < 1000) {
                    print_time_and_process(next_time,running->name);
                    std::cout << " (tau " << running->tau << "ms) completed a CPU burst; " << (running->cpu_burst_time).size() << " burst";
                    if ((running->cpu_burst_time).size() != 1) std::cout << "s";
                    std::cout << " to go";
                    print_priority_queue(ready_queue);
                }

                int new_tau = recalculate_tau(running->tau,runtime,alpha);

                if (next_time < 1000) {
                    std::cout << "time " << next_time << "ms: Recalculated tau for process " << running->name << ": old tau " << running->tau << "ms; new tau " << new_tau << "ms";
                    print_priority_queue(ready_queue);
                    print_time_and_process(next_time,running->name);
                    std::cout << " switching out of CPU; will block on I/O until time " << next_time + context_time/2 + (running->io_burst_time).front() << "ms";
                    print_priority_queue(ready_queue);
                }
                running->tau = new_tau;
                running->next_event_time = next_time + context_time/2 + (running->io_burst_time).front();
                io_wait.push(*running);

            }
            else {
                print_time_and_process(next_time,running->name);
                std::cout << " terminated";
                print_priority_queue(ready_queue);
            }

            delete(running);
            running = nullptr;
            cpu_ready_time = next_time + context_time;
            choose_next_process = next_time + context_time/2;
        }

        else if (next_event == "entering cpu burst") {
            context_switches++;
            wait_time += next_time - context_time/2;
            if (next_time < 1000) {
                print_time_and_process(next_time,running->name);
                std::cout << " (tau " << running->tau << "ms) started using the CPU for ";
                if (running->cpu_time_left != (running->cpu_burst_time).front()) {
                    std::cout << "remaining " << running->cpu_time_left << "ms of ";
                }
                std::cout << (running->cpu_burst_time).front() << "ms burst";
            }
            if (next_time < 1000) {
                print_priority_queue(ready_queue);
            }
            if (ready_queue.empty()) {
                running->next_event_time = next_time + running->cpu_time_left;
                choose_next_process = INT_MAX;
            }
            else {
                Process p = ready_queue.top();
                if (running->tau - ((running->cpu_burst_time).front() - running->cpu_time_left) >
                (p.tau - (p.cpu_burst_time.front() - p.cpu_time_left))) {
                    if (next_time < 1000) {
                        print_time_and_process(next_time,p.name);
                        std::cout << " (tau " << p.tau << "ms) will preempt " << running->name;
                        print_priority_queue(ready_queue);
                    }
                    preemptions++;
                    wait_time -= next_time + context_time/2;
                    choose_next_process = next_time + context_time/2;
                    cpu_ready_time = next_time + context_time;
                    ready_queue.push(*running);
                    delete(running);
                    running = nullptr;
                }
                else {
                    running->next_event_time = next_time + running->cpu_time_left;
                    choose_next_process = INT_MAX;
                }
            }
        }

        else if (next_event == "io burst complete") {
            Process p = io_wait.top();
            io_wait.pop();
            p.io_burst_time.pop_front();
            if (ready_queue.empty() && running == nullptr && cpu_ready_time < next_time + context_time/2) {
                choose_next_process = next_time;
                cpu_ready_time = next_time + context_time/2;
            }
            p.ready_time = next_time;
            wait_time -= next_time;
            turnaround_time -= next_time;
            ready_queue.push(p);

            if (next_time < 1000) {
                print_time_and_process(next_time,p.name);
                std::cout << " (tau " << p.tau << "ms) completed I/O; ";
            }
            if (running != nullptr && p.tau < (running->tau - ((running->cpu_burst_time).front() - (running->next_event_time - next_time)))) {
                if (next_time < 1000) {
                    std::cout << "preempting " << running->name;
                    print_priority_queue(ready_queue);
                }
                preemptions++;
                wait_time -= next_time + context_time/2;
                turnaround_time += next_time - running->ready_time + context_time / 2;
                running->cpu_time_left = running->next_event_time - next_time;
                ready_queue.push(*running);
                delete(running);
                running = nullptr;
                if (choose_next_process == INT_MAX) {
                    choose_next_process = next_time + context_time/2;
                    cpu_ready_time = next_time + context_time;
                }
            }
            else {
                if (next_time < 1000) {
                    std::cout << "added to ready queue";
                    print_priority_queue(ready_queue);
                }
            }
        }

        else if (next_event == "process enters the scheduler") {
            Process p = unstarted.top();
            p.cpu_time_left = p.cpu_burst_time.front();
            unstarted.pop();
            if (ready_queue.empty() && running == nullptr && cpu_ready_time < next_time + context_time/2) {
                choose_next_process = next_time;
                cpu_ready_time = next_time + context_time/2;
            }
            wait_time -= next_time;
            turnaround_time -= next_time;
            ready_queue.push(p);
            if (next_time < 1000) {
                print_time_and_process(next_time,p.name);
                std::cout << " (tau " << p.tau << "ms) arrived; added to ready queue";
                print_priority_queue(ready_queue);
            }
        }
        else {
            std::cout << "Unknown event occurred." << std::endl;
        }
    }
    total_time = next_time + context_time / 2;
    std::cout << "time " << next_time + context_time/2 << "ms: Simulator ended for SRT [Q: empty]" << std::endl;
    print_algorithm_stats(output,total_number_bursts,wait_time,turnaround_time,total_bursts_time,
                            context_switches,preemptions,total_time,0,0,1,0);
}

void rr(const std::vector<Process>& process_list, int context_time, int preemption_time, std::ofstream& output) {
    Process* running = nullptr;
    std::list<Process> ready_queue;
    std::priority_queue<Process> io_wait;
    std::priority_queue<Process> unstarted;
    int cpu_ready_time = context_time/2;
    int next_time = 0;
    int finish_cpu_time = INT_MAX;
    int choose_next_process = 0;

    int total_bursts_time = 0;
    int total_number_bursts = 0;
    int wait_time = 0;
    //int arrival_time;
    int turnaround_time = 0;
    int context_switches = 0;
    int preemptions = 0;
    int total_time = 0;
    Process* next_process = nullptr;
    for (unsigned int i = 0; i < process_list.size(); i++) {
        unstarted.push(process_list[i]);
        total_number_bursts += process_list[i].cpu_burst_time.size();
    }
    std::cout << "time 0ms: Simulator started for RR with time slice " << preemption_time << "ms [Q: empty]" << std::endl;
    while (running != nullptr || ready_queue.size() > 0 || io_wait.size() > 0 || unstarted.size() > 0) {
        next_time = INT_MAX;
        std::string next_event = "";
        if (running != nullptr) {
            next_time = running->next_event_time;
            next_event = "cpu burst complete";
        }
        if ((running == nullptr || finish_cpu_time != INT_MAX) && !ready_queue.empty()) {
            if (cpu_ready_time < next_time) {
                next_time = cpu_ready_time;
                next_event = "entering cpu burst";
            }
        }
        if (!io_wait.empty()) {
            Process p = io_wait.top();
            if (p.next_event_time < next_time) {
                next_time = p.next_event_time;
                next_event = "io burst complete";
            }
        }
        if (!unstarted.empty()) {
            Process p = unstarted.top();
            if (p.next_event_time < next_time) {
                next_time = p.next_event_time;
                next_event = "process enters the scheduler";
            }
        }

        if (next_time >= finish_cpu_time) {
            ready_queue.push_back(*running);
            delete(running);
            running = nullptr;
            finish_cpu_time = INT_MAX;
        }
        if (next_time > choose_next_process && running == nullptr && !ready_queue.empty() && next_process == nullptr) {
            next_process = new Process();
            *next_process = ready_queue.front();
            next_process->next_event_time = choose_next_process + context_time/2 + std::min(preemption_time,next_process->cpu_time_left);
            ready_queue.pop_front();
        }

        if (next_event == "cpu burst complete") {
            if (running->cpu_time_left > preemption_time) {
                running->cpu_time_left -= preemption_time;
                if (ready_queue.empty()) {
                    if (next_time < 1000) {
                        std::cout << "time " << next_time << "ms: Time slice expired; no preemption because ready queue is empty";
                        print_list_queue(ready_queue);
                    }
                    running->next_event_time = next_time + std::min(preemption_time,running->cpu_time_left);
                }
                else {
                    preemptions++;
                    wait_time -= next_time + context_time/2;
                    if (next_time < 1000) {
                        std::cout << "time " << next_time << "ms: Time slice expired; process " << running->name << " preempted with " << running->cpu_time_left << "ms remaining";
                        print_list_queue(ready_queue);
                    }
                    finish_cpu_time = next_time + context_time/2;
                    choose_next_process = next_time + context_time/2;
                    cpu_ready_time = next_time + context_time;
                    running->next_event_time = INT_MAX;
                }
            }
            else {
                total_bursts_time += running->cpu_burst_time.front();
                (running->cpu_burst_time).pop_front();

                if ((running->io_burst_time).size() > 0) {
                    running->cpu_time_left = (running->cpu_burst_time).front();
                    if (next_time < 1000) {
                        print_time_and_process(next_time,running->name);
                        std::cout << " completed a CPU burst; " << (running->cpu_burst_time).size() << " burst";
                        if ((running->cpu_burst_time).size() != 1) std::cout << "s";
                        std::cout << " to go";
                        print_list_queue(ready_queue);
                        print_time_and_process(next_time,running->name);
                        std::cout << " switching out of CPU; will block on I/O until time " << next_time + context_time/2 + (running->io_burst_time).front() << "ms";
                        print_list_queue(ready_queue);
                    }

                    turnaround_time += next_time - running->ready_time + context_time / 2;
                    running->next_event_time = next_time + context_time/2 + (running->io_burst_time).front();
                    io_wait.push(*running);
                }
                else {
                    turnaround_time += next_time - running->ready_time + context_time / 2;
                    print_time_and_process(next_time,running->name);
                    std::cout << " terminated";
                    print_list_queue(ready_queue);
                }

                delete(running);
                running = nullptr;
                choose_next_process = next_time + context_time/2;
                cpu_ready_time = next_time + context_time;
            }
        }

        else if (next_event == "entering cpu burst") {
            running = next_process;
            next_process = nullptr;
            wait_time += next_time - context_time / 2;
            context_switches++;
            if (next_time < 1000) {
                print_time_and_process(next_time,running->name);
                std::cout << " started using the CPU for ";
                if (running->cpu_time_left != (running->cpu_burst_time).front()) {
                    std::cout << "remaining " << running->cpu_time_left << "ms of ";
                }
                std::cout << (running->cpu_burst_time).front() << "ms burst";
                print_list_queue(ready_queue);
            }
        }

        else if (next_event == "io burst complete") {
            Process p = io_wait.top();
            io_wait.pop();
            p.io_burst_time.pop_front();
            if (ready_queue.empty() && cpu_ready_time < next_time + context_time/2 && next_process == nullptr) {
                choose_next_process = next_time;
                cpu_ready_time = next_time + context_time/2;
            }
            p.ready_time = next_time;
            wait_time -= next_time;
            ready_queue.push_back(p);

            if (next_time < 1000) {
                print_time_and_process(next_time,p.name);
                std::cout << " completed I/O; added to ready queue";
                print_list_queue(ready_queue);
            }
        }

        else if (next_event == "process enters the scheduler") {
            Process p = unstarted.top();
            p.cpu_time_left = p.cpu_burst_time.front();
            unstarted.pop();
            if (ready_queue.empty() && cpu_ready_time < next_time + context_time/2 && next_process == nullptr) {
                choose_next_process = next_time;
                cpu_ready_time = next_time + context_time/2;
            }
            wait_time -= next_time;
            p.ready_time = next_time;
            ready_queue.push_back(p);

            if (next_time < 1000) {
                print_time_and_process(next_time,p.name);
                std::cout << " arrived; added to ready queue";
                print_list_queue(ready_queue);
            }
        }
        else {
            std::cout << "Unknown event occurred." << std::endl;
        }
    }
    total_time = next_time + context_time/2;
    std::cout << "time " << next_time + context_time/2 << "ms: Simulator ended for RR [Q: empty]" << std::endl;
    print_algorithm_stats(output,total_number_bursts,wait_time,turnaround_time,total_bursts_time,
                            context_switches,preemptions,total_time,0,0,0,1);
}



int main(int argc, char* argv[]) {

	if(argc != 8)
	{
		std::cerr << "ERROR: Incorrect Number of Arguments" << std::endl;
		return EXIT_FAILURE;
	}

	int n = atoi(argv[1]);

	if(n < 1 || n > 26)
	{
		std::cerr << "ERROR: Invalid Number of Processes" << std::endl;
		return EXIT_FAILURE;
	}

	if(!isdigit(argv[2][0]) && !(argv[2][0] == '-' && isdigit(argv[2][1])))
	{
			std::cerr << "ERROR: Invalid Seed" << std::endl;
			return EXIT_FAILURE;
	}

	int seed = atoi(argv[2]);

	double lambda = atof(argv[3]);

	if(lambda <= 0.0)
	{
		std::cerr << "ERROR: Invalid Lambda Value" << std::endl;
		return EXIT_FAILURE;
	}

	double upperBound = atof(argv[4]);

	if(upperBound <= 0)
	{
		std::cerr << "ERROR: Invalid Upper Bound" << std::endl;
		return EXIT_FAILURE;
	}

	int Tcs = atoi(argv[5]);

	if(Tcs <= 0 || Tcs % 2)
	{
		std::cerr << "ERROR: Invalid Context Switch Time" << std::endl;
		return EXIT_FAILURE;
	}

	if(!isdigit(argv[6][0]) && !(argv[6][0] == '-' && isdigit(argv[6][1])))
	{
			std::cerr << "ERROR: Alpha Value" << std::endl;
			return EXIT_FAILURE;
	}

	double alpha = atof(argv[6]);

	int Ts = atoi(argv[7]);

	if(Ts <= 0)
	{
		std::cerr << "ERROR: Invalid Slice Time" << std::endl;
		return EXIT_FAILURE;
	}

	std::vector<Process> lst;

	srand48(seed);

	for(char i = 'A'; i < 'A' + n; i++)
	{
		Process p;
		p.name = i;
		p.tau = 1/lambda;

		p.next_event_time = floor(next_exp(lambda, upperBound));

		int bursts = ceil(drand48()*100);

		std::cout << "Process " << i << ": arrival time " << p.next_event_time << "ms; tau " << p.tau << "ms; " << bursts << " CPU bursts:" << std::endl;

		for(int k = 1; k < bursts; k++)
		{
			int cpu = ceil(next_exp(lambda, upperBound));
			int io = ceil(next_exp(lambda, upperBound))*10;

			std::cout << "--> CPU burst " << cpu << "ms --> I/O burst " << io << "ms" << std::endl;

			p.cpu_burst_time.push_back(cpu);
			p.io_burst_time.push_back(io);
		}

		int cpu = ceil(next_exp(lambda, upperBound));
		std::cout << "--> CPU burst " << cpu << "ms" << std::endl;
		p.cpu_burst_time.push_back(cpu);

		lst.push_back(p);
	}

	std::cout << std::endl;

	std::ofstream statistics("simout.txt");
    fcfs(lst,Tcs,statistics);
    std::cout << "\n";
    sjf(lst,Tcs,alpha,statistics);
    std::cout << "\n";
    srt(lst,Tcs,alpha,statistics);
    std::cout << "\n";
    rr(lst,Tcs,Ts,statistics); //*/

}
