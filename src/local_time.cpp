#include <sys/stat.h>
#include <sys/time.h>
#include <stdio.h>
class LocalTime {
public:
    unsigned short year, month, date, day_of_week;
    unsigned short hour, minute, second, millisecond;
public:
    void get(){
        time_t ct; ct=time(&ct); struct tm * st = (struct tm *)localtime(&ct);
        year    = st->tm_year + 1900;
        month   = st->tm_mon;
        date    = st->tm_mday;
        day_of_week = st->tm_wday;
        hour    = st->tm_hour;
        minute  = st->tm_min;
        second  = st->tm_sec;
        millisecond = 0;
    }
    const char * month_name(){
        switch (month){
            case 0 : return "Jan"; case 1 : return "Feb";
            case 2 : return "Mar"; case 3 : return "Apr";
            case 4 : return "May"; case 5 : return "Jun";
            case 6 : return "Jul"; case 7 : return "Aug";
            case 8 : return "Sep"; case 9 : return "Oct";
            case 10: return "Nov"; case 11: return "Dec";
            default:    return "???";
        }
    }
    void translate(char * buffer){ sprintf(buffer, "%s %d, %4d %02d:%02d", month_name(), date, year, hour, minute); }
    void translate_precise(char * buffer){ sprintf(buffer, "%4d.%02d.%02d %02d:%02d:%02d", year, month, date, hour, minute, second); }
    void translate_ms_precise(char * buffer){ sprintf(buffer, "%4d.%02d.%02d %02d:%02d:%02d.%03d", year, month, date, hour, minute, second, millisecond); }
    void get_text(char * buffer){ get(); translate(buffer); }
    void get_text_precise(char * buffer){ get(); translate_precise(buffer); }
    void get_text_ms_precise(char * buffer){ get(); translate_ms_precise(buffer); }
};
