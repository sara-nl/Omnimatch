void tick(struct tm *start)
{
time_t lt;
lt=time(NULL);
start=localtime(&lt);
/* end tick */
}
void tack(struct tm *start)
{
struct tm *stop;
time_t lt;
lt=time(NULL);
stop=localtime(&lt);
printf("Time : %i:%i:%i\n",stop->tm_hour,stop->tm_min,stop->tm_sec);

/* end tack */
}
