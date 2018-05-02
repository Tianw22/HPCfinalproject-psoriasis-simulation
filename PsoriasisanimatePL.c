# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <math.h>
# include <time.h>
# include <sys/time.h>
# include <unistd.h>
# include <omp.h>

# define NX 40
# define NY 40
#define ANSI_COLOR_RED     "\x1b[41m"
#define ANSI_COLOR_GREEN   "\x1b[42m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_LTGRAY   "\x1b[100m"
#define ANSI_COLOR_RESET   "\x1b[0m"
#define CLEAR_SCREEN       "\e[1;1H\e[2J"

typedef struct {
    char STATE;
    double B,I,D;
    int i,j;
} psoriasis;

void timestamp ( void );
void showpsoriasis (int, int, psoriasis**);
void showpsoriasis_persist (int, int, psoriasis**);
void p(const char *str);


int main ( void )
{
    int i;
    int j;
    int nx = NX;
    int ny = NY;
    int NSTEPS = 100;
    struct timeval start_time, stop_time, elapsed_time;

    psoriasis **t,**tnew;

    t = malloc(nx*sizeof(psoriasis *));
    *t = malloc(nx*ny*sizeof(psoriasis));
    for (i = 0; i < nx; i++)
        *(t+i) = *t + i*ny;

    tnew = malloc(nx*sizeof(psoriasis *));
    *tnew = malloc(nx*ny*sizeof(psoriasis));
    for (i = 0; i < nx; i++)
        *(tnew+i) = *tnew + i*ny;

    srand(time(NULL));
    double Dinit=1;

    gettimeofday(&start_time,NULL);

    omp_set_num_threads(16);
    /* Initialize skin -- no lesion yet */
    for ( j = 0; j < ny; j++ )
    {
        for ( i = 0; i < nx; i++ )
        {
            // Cell sites are aware of their own location.
            tnew[i][j].i = i;
            tnew[i][j].j = j;

            // No lesion around edge of skin
            if( i == 0 || j == 0 || i == nx -1 || j == ny-1 )
            {
                tnew[i][j].STATE = ' ';
            } else {
                // populate skin with psoriasis
                t[i][j].D = tnew[i][j].D = (double) rand() /RAND_MAX;
                t[i][j].B = tnew[i][j].B = 0.5;// (double) rand() /RAND_MAX;
                t[i][j].I = tnew[i][j].I = 0.95;//(double) rand() /RAND_MAX;
                if ( Dinit - tnew[i][j].D > 0 ) {
                    tnew[i][j].STATE = '^';
                } else {
                    tnew[i][j].STATE = ' ';
                }
            }
        }
    }

    /* Start a lesion in the middle of the grid */

    /* Constant source lesion */ 
    tnew[nx/2][ny/2].STATE = 'D';

    showpsoriasis(nx,ny,tnew);

    /* Let it develop */

    int it;

    for (it=0 ;it<NSTEPS ;it++ )
    {

        /* Save the current lesion state. */
      
        for ( j = 0; j < ny; j++ )
        {
            #pragma omp parallel for private(i) reduction(max:it) schedule(dynamic)
            for ( i = 0; i < nx; i++ )
            {
                t[i][j].STATE = tnew[i][j].STATE;
            }
        }

        /*Scan for damaged neighbors.*/

        for ( j = 1; j < ny-1; j++ )
        {
            #pragma omp parallel for private(i) reduction(max:it) schedule(dynamic)
            for ( i = 1; i < nx-1; i++ )
            {
                //If a cell is damaged see if it continues developing
                //otherwise the lesion goes out.
               
                if ( t[i][j].STATE == 'D' )
                {
                    if (t[i][j].B < (double) rand() /RAND_MAX)
                    {
                        tnew[i][j].STATE = '.';
                    }
                }

                // If cell is health but has damaged neighbors see if cell gets lesion.
                if ( t[i][j].STATE == '^' )
                {
                    
                    //It will stop growing when it develop into a specific size.
                   if (((20-i)*(20-i)+(20-j)*(20-j))<225)
                   {
                    // Either a corner neighbor possibly infects cell or side neighbor;
                    // Corner neighbor influence is suppressed but 1/sqrt(2) 0.293 is approx = 1-1/sqrt(2)
                    if ( 1*0.293 > (double) rand()/ RAND_MAX)
                    {
                        if (t[i-1][j-1].STATE == 'D' ||
                                t[i-1][j+1].STATE == 'D' ||
                                t[i+1][j-1].STATE == 'D' ||
                                t[i+1][j+1].STATE == 'D' )
                        {
                            if (t[i][j].I > (double) rand() /RAND_MAX)
                            {
                                tnew[i][j].STATE = 'D';
                            }
                        }
                    } else {
                        if (t[i-1][j].STATE == 'D' ||
                                t[i][j-1].STATE == 'D' ||
                                t[i][j+1].STATE == 'D' ||
                                t[i+1][j].STATE == 'D' )
                        {
                            if (t[i][j].I > (double) rand() /RAND_MAX)
                            {
                                tnew[i][j].STATE = 'D';
                            }
                        }
                    }
                   }

                }

                //Recover from middle.

                if ( t[i][j].STATE == '.' )
                {
                    if (
                        t[(i-2) > 0 ? i-2 : 0][(j-2) > 0 ? j-2 : 0].STATE != 'D' &&
                        t[(i-2) > 0 ? i-2 : 0][(j-1) > 0 ? j-1 : 0].STATE != 'D' &&
                        t[(i-2) > 0 ? i-2 : 0][j].STATE != 'D' &&
                        t[(i-2) > 0 ? i-2 : 0][(j+1) > ny-1 ? ny-1 : j+1].STATE != 'D' &&
                        t[(i-2) > 0 ? i-2 : 0][(j+2) > ny-1 ? ny-1 : j+2].STATE != 'D' &&

                        t[(i-1) > 0 ? i-1 : 0][(j-2) > 0 ? j-2 : 0].STATE != 'D' &&
                        t[(i-1) > 0 ? i-1 : 0][(j-1) > 0 ? j-1 : 0].STATE != 'D' &&
                        t[(i-1) > 0 ? i-1 : 0][j].STATE != 'D' &&
                        t[(i-1) > 0 ? i-1 : 0][(j+1) > ny-1 ? ny-1 : j+1].STATE != 'D' &&
                        t[(i-1) > 0 ? i-1 : 0][(j+2) > ny-1 ? ny-1 : j+2].STATE != 'D' &&

                        t[i][(j-2) > 0 ? j-2 : 0].STATE != 'D' &&
                        t[i][(j-1) > 0 ? j-1 : 0].STATE != 'D' &&
                        t[i][(j+1) > ny-1 ? ny-1 : j+1].STATE != 'D' &&
                        t[i][(j+2) > ny-1 ? ny-1 : j+2].STATE != 'D' &&

                        t[(i+1) > nx-1 ? nx-1 : i+1][(j-2) > 0 ? j-2 : 0].STATE != 'D' &&
                        t[(i+1) > nx-1 ? nx-1 : i+1][(j-1) > 0 ? j-1 : 0].STATE != 'D' &&
                        t[(i+1) > nx-1 ? nx-1 : i+1][j].STATE != 'D' &&
                        t[(i+1) > nx-1 ? nx-1 : i+1][(j+1) > ny-1 ? ny-1 : j+1].STATE != 'D' &&
                        t[(i+1) > nx-1 ? nx-1 : i+1][(j+2) > ny-1 ? ny-1 : j+2].STATE != 'D' &&

                        t[(i+2) > nx-1 ? nx-1 : i+2][(j-2) > 0 ? j-2 : 0].STATE != 'D' &&
                        t[(i+2) > nx-1 ? nx-1 : i+2][(j-1) > 0 ? j-1 : 0].STATE != 'D' &&
                        t[(i+2) > nx-1 ? nx-1 : i+2][j].STATE != 'D' &&
                        t[(i+2) > nx-1 ? nx-1 : i+2][(j+1) > ny-1 ? ny-1 : j+1].STATE != 'D' &&
                        t[(i+2) > nx-1 ? nx-1 : i+2][(j+2) > ny-1 ? ny-1 : j+2].STATE != 'D'
                    )
                    {
                        tnew[i][j].STATE = 'o';
                    }

                }

            }

        }
        showpsoriasis(nx,ny,tnew);
    }
    gettimeofday(&stop_time,NULL);

    puts("AFTER RECOVERY");
    printf ( "\n" );
    printf ( "Almeida et al.,\n \tJournal of Physics: \n\tConference Series 285 (2011) 012038 \n\tdoi:10.1088/1742-6596/285/1/012038:\n" );
    printf ("Modeling Pattern Formation in Skin Diseases by a Cellular Automaton \n \t Journal of Investigative Dermatology (2013) 133, 567–571\n\t doi:10.1038/jid.2012.321 \n\t published online 30 August 2012");
    printf ( "\n" );
    timestamp ( );
    timersub(&stop_time, &start_time, &elapsed_time);

    free(*t);
    free(*tnew);
    free(t);
    free(tnew);
    return 0;
}
/******************************************************************************/

void timestamp ( void )
{
# define TIME_SIZE 40

    static char time_buffer[TIME_SIZE];
    const struct tm *tm;
    time_t now;

    now = time ( NULL );
    tm = localtime ( &now );

    strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

    printf ( "%s\n", time_buffer );

    return;
# undef TIME_SIZE
}

/******************************************************************************/
void showpsoriasis(int nx, int ny, psoriasis **tnew)
{
    int i,j;
    for ( j = 0; j < ny; j++ )
    {
        for ( i = 0; i < nx; i++ )
        {
            if (tnew[i][j].STATE == 'D')
            {
                printf(ANSI_COLOR_MAGENTA"%c "ANSI_COLOR_RESET,tnew[i][j].STATE);
            }
            else if (tnew[i][j].STATE == '^')
            {
                printf(ANSI_COLOR_GREEN"%c "ANSI_COLOR_RESET,tnew[i][j].STATE);
            }else if (tnew[i][j].STATE == '.')
            {
                printf(ANSI_COLOR_RED"%c "ANSI_COLOR_RESET,tnew[i][j].STATE);
            }else if (tnew[i][j].STATE == 'o')
            {
                printf(ANSI_COLOR_GREEN"%c "ANSI_COLOR_RESET,tnew[i][j].STATE);

            } else {
                printf("%c ",tnew[i][j].STATE);
            }


        }
        printf("\n");
    }

    sleep(1);
    printf(CLEAR_SCREEN);

}
/******************************************************************************/
void showpsoriasis_persist(int nx, int ny, psoriasis **tnew)
{
    int i,j;
    for ( j = 0; j < ny; j++ )
    {
        for ( i = 0; i < nx; i++ )
        {
            if (tnew[i][j].STATE == 'D')
            {
                printf(ANSI_COLOR_MAGENTA"%c "ANSI_COLOR_RESET,tnew[i][j].STATE);
            }
            else if (tnew[i][j].STATE == '^')
            {
                printf(ANSI_COLOR_GREEN"%c "ANSI_COLOR_RESET,tnew[i][j].STATE);
            }else if (tnew[i][j].STATE == '.')
            {
                printf(ANSI_COLOR_RED"%c "ANSI_COLOR_RESET,tnew[i][j].STATE);
            }else if (tnew[i][j].STATE == 'o')
            {
                printf(ANSI_COLOR_GREEN"%c "ANSI_COLOR_RESET,tnew[i][j].STATE);

            } else {
                printf("%c ",tnew[i][j].STATE);
            }


        }
        printf("\n");
    }

}

/******************************************************************************/
