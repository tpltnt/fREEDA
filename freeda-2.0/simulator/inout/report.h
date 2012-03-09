/**************************************************************************
 *
 * Use this function to send output messages.
 *
 **************************************************************************/

#ifndef report_h
#define report_h

/* Usage:
 * LOGFILE puts 'message' in log file
 * MESSAGE puts 'message' on screen display and in log file
 * WARNING puts 'message' on screen display and in log file
 * FATAL puts 'message' on screen display and in log file, exits */

enum ReportType {LOGFILE, MESSAGE, WARNING, FATAL};

void report(enum ReportType type, const char *message);

void sepLine();

void newLine();

#endif
