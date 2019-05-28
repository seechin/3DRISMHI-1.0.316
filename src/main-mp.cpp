#ifdef _LOCALPARALLEL_
    void subroutine_main(IET_Param * sys, IET_arrays * arr, int id){
        pid_t ppid = sys->pid; bool mp_by_fork = sys->mp_by_fork;
        pid_t pidthis = __forkpid?__forkpid[id]:ppid;
        #ifdef _LOCALPARALLEL_PTHREAD_
            __thread_list[id] = pthread_self();
        #else
            __thread_list[id] = 0;
        #endif
        int jobs_done_total = 0; double jobs_time_total = 0;
        //double get_current_time_double

        if (id<=0) return;
        int ndeadclk = 0;
        bool alert_on_exit = sys->debug_level>=1? true : false;
        //fprintf(stderr, "FORK: %d\n", id); fflush(stderr);
        while (true){ double time_stamp = get_current_time_double();
            if (sys->mp_tasks[id] == MPTASK_TERMINATE){
                sys->mp_tasks[id] = MPTASK_NONE;
                break;
            } else if (sys->mp_tasks[id] == MPTASK_NONE){
                usleep(100);
                if (mp_by_fork){
                    ndeadclk ++; if (ndeadclk > 1000){
                        kill(ppid, 0); if (errno == ESRCH){
                            if (alert_on_exit) fprintf(stderr, "%s : error : pid %d will end as parent pid %d disappeared\n", software_name, pidthis, ppid);
                            alert_on_exit = false; break;
                        }; ndeadclk = 0;
                    }
                }
                continue;
            } else if (sys->mp_tasks[id] == MPTASK_FFTW){
                perform_3rx1k_convolution_r1(id, &arr->fftw_mp, arr->fftw_mp.fft[id].si, arr->fftw_mp.fft[id].sj);
                jobs_done_total ++; jobs_time_total += get_current_time_double() - time_stamp;
            } else if (sys->mp_tasks[id] == MPTASK_MERGE_FFT_DATA){
                sumback_3rx1k_convolution_t1(&arr->fftw_mp, id, arr->fftw_mp.n_active_jobs);
            } else if (sys->mp_tasks[id] == MPTASK_HI_SOLVER){
                subroutine_hi_solver(sys, arr, id);
                jobs_done_total ++; jobs_time_total += get_current_time_double() - time_stamp;
            } else if (sys->mp_tasks[id] == MPTASK_DIIS_STEPIN){
                arr->diis_current->step_in_mp_parallel(id);
                jobs_done_total ++; jobs_time_total += get_current_time_double() - time_stamp;
            } else if (sys->mp_tasks[id] == MPTASK_DIIS_WEIGHT){
                arr->diis_current->calc_weights_sub(id);
                jobs_done_total ++; jobs_time_total += get_current_time_double() - time_stamp;
            } else if (sys->mp_tasks[id] == MPTASK_RISM_CLOSURE){
                RISMHI3D_RISMNS::perform_closure(sys, arr, id);
                jobs_done_total ++; jobs_time_total += get_current_time_double() - time_stamp;
            } else if (sys->mp_tasks[id] == MPTASK_FFSR){
                build_force_field_sr_1(arr->ffsr_mp.irange[id][0], arr->ffsr_mp.irange[id][1], &arr->ffsr_mp.irange[id][2], sys, arr, arr->ffsr_mp.lj[id], arr->ffsr_mp.coulsr[id], arr->ffsr_mp.r2uvmin[id], arr->ffsr_mp.coulp2[id], arr->ffsr_mp.pseudoliquid_potential[id]);
                jobs_done_total ++; jobs_time_total += get_current_time_double() - time_stamp;
            } else if (sys->mp_tasks[id] == MPTASK_MERGE_FF_DATA){
                merge_force_field_mp_data(sys, arr, id);
                jobs_done_total ++; jobs_time_total += get_current_time_double() - time_stamp;
            }
            sys->mp_tasks[id] = MPTASK_NONE;
        }
        join_here(id, mp_by_fork);
        if (alert_on_exit && sys->log()){ char buffer[1024];
            if (id==0) fprintf(sys->log(), "DEBUG:: %s ends normally\n", mp_by_fork?"process":"thread[0]");
            else if (sys->log()==stdout||sys->log()==stderr) fprintf(sys->log(), "DEBUG:: %s[%d] done after %d jobs, %s sec\n", mp_by_fork?"process":"thread", id, jobs_done_total, display_time(jobs_time_total / 1000.0, buffer));
        }
        #ifdef _LOCALPARALLEL_PTHREAD_
            pthread_exit(nullptr);
        #else
            exit(0);
        #endif
    }

    // Following functions are called only by main thread
    bool wait_subroutines(int nt, int mp_tasks[MAX_THREADS], int timeup_ms){
        if (nt<=1) return true;
        double time0 = get_current_time_double();
        while (true){
            bool finished = true; for (int i=1; i<nt; i++) if (mp_tasks[i]!=MPTASK_NONE){ finished = false; break; }
            if (finished) return true;
            if (timeup_ms>0 && get_current_time_double()-time0 > timeup_ms) return false;
            usleep(100);
        }
        return true;
    }
    bool wait_subroutines(IET_Param * sys, int timeup_ms){
        return wait_subroutines(sys->nt, sys->mp_tasks, timeup_ms);
    }
    bool wait_subroutines_end(IET_Param * sys, int timeup_ms=-1){
        if (sys->nt<=1) return true;
        /*while (true){
            char buffer[MAX_THREADS+10]; memset(buffer, 0, sizeof(buffer)); buffer[0] = ':';
            for (int i=1; i<sys->nt; i++){
                //buffer[i] = sys->mp_tasks[i]==MPTASK_NONE? '.' : '+';
                if (sys->mp_by_fork){
                    kill(__forkpid[i], 0); if (errno == ESRCH) buffer[i] = '.'; else buffer[i] = '+';
                } else {
                    pthread_kill(__thread_list[i], 0); if (errno == ESRCH) buffer[i] = '.'; else buffer[i] = '+';
                }
                if (buffer[i]=='+' && sys->mp_tasks[i]==MPTASK_NONE) buffer[i] = '?';
            }
            static int iii; iii ++;
            fprintf(stderr, "%3d %s\r", iii%1000, buffer);
            usleep(1000);
        }*/
        double time0 = get_current_time_double();
        while (true){
            int alive = 0; for (int i=1; i<sys->nt; i++){
              #ifdef _LOCALPARALLEL_PTHREAD_
                if (sys->mp_by_fork){
                    if (sys->mp_tasks[i]==MPTASK_NONE){
                    } else {
                        kill(__forkpid[i], 0);
                        if (errno != ESRCH) alive ++;
                    }
                } else {
                    pthread_kill(__thread_list[i], 0);
                    if (errno != ESRCH) alive ++;
                }
              #else
                if (sys->mp_tasks[i]==MPTASK_NONE){
                } else {
                    kill(__forkpid[i], 0);
                    if (errno != ESRCH) alive ++;
                }
              #endif
            }
            if (alive<1) return true;
            if (timeup_ms>0 && get_current_time_double()-time0 > timeup_ms) return false;
            usleep(100);
        }
        return true;
    }

  #ifdef _LOCALPARALLEL_PTHREAD_
    class pThread_subroutine_param { public: int id; IET_Param * sys; IET_arrays * arr; };
    pThread_subroutine_param __thread_params[MAX_THREADS];
    void * thread_entry(void * ptp){
        pThread_subroutine_param * pt = (pThread_subroutine_param*) ptp;
        //printf(" thread %d of %d\n", pt->id, pt->sys->nt); return nullptr;
        subroutine_main(pt->sys, pt->arr, pt->id);
        return nullptr;
    }
    void create_subroutines_by_pthread(IET_Param * sys, IET_arrays * arr){
        for (int i=0; i<MAX_THREADS; i++){ __thread_params[i].id = i; __thread_params[i].sys = sys; __thread_params[i].arr = arr; }
        if (sys->nt<=1) return;
        for (int i=1; i<sys->nt; i++){
            int ptret = pthread_create(&__thread_list[i], nullptr, thread_entry, &__thread_params[i]);
        }
    }
    void kill_subroutines_by_pthread(IET_Param * sys){
        for (int i=1; i<sys->nt; i++){
            sys->mp_tasks[i] = MPTASK_TERMINATE;
            pthread_kill(__thread_list[i], SIGKILL/*9*/);
        }
    }
  #endif

    void create_subroutines_by_fork(IET_Param * sys, IET_arrays * arr){
        if (sys->nt<=1) return;
        int id = fork_here(sys->nt);
        subroutine_main(sys, arr, id);
    }
    void kill_subroutines_by_fork(IET_Param * sys){
        for (int i=1; i<sys->nt; i++){
            sys->mp_tasks[i] = MPTASK_TERMINATE;
            kill(__forkpid[i], SIGKILL/*9*/);
        }
    }

    void create_subroutines(IET_Param * sys, IET_arrays * arr){
      #ifdef _LOCALPARALLEL_PTHREAD_
        if (sys->mp_by_fork){
            create_subroutines_by_fork(sys, arr);
        } else {
            create_subroutines_by_pthread(sys, arr);
        }
      #else
        create_subroutines_by_fork(sys, arr);
      #endif
    }
    void kill_subroutines(IET_Param * sys){
      #ifdef _LOCALPARALLEL_PTHREAD_
        if (sys->mp_by_fork){
            kill_subroutines_by_fork(sys);
        } else {
            kill_subroutines_by_pthread(sys);
        }
      #else
        kill_subroutines_by_fork(sys);
      #endif
    }

#endif
