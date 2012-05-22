/*
 *
 * searching of a (n, (n+1)/2) code with the given distance d.
 *
 */

#include <stdio.h>
#include <sys/vtimes.h>

typedef unsigned short uShort;

#define exp2(n) (1<<(n))
#define mask(n) (exp2(n))
#define weight(c) (dependTab[0][c])
#define min(x,y) ((x)<(y)?(x):(y))
#define max(x,y) ((x)>(y)?(x):(y))

#define MAXK 16
#define TRUE 1
#define FALSE 0
#define CRITERIA 0

char dependTab[MAXK][exp2(MAXK)];

uShort pool[MAXK][exp2(MAXK)];
uShort *poolptr[MAXK + 1];
uShort rowValue[MAXK + 1][MAXK + 1];

int n, k, m, d;
int poolnum[MAXK];

struct vtimes vtms;

extern void matrixInit(), printVec(), printResult(), calWeight();

void main(argc, argv)
  int argc;
  char *argv[];
{
  int ac = argc - 1;
  char **av = argv + 1;
 
  n = atoi(av[0]);
  d = atoi(av[1]);
  if (n % 2 != 1) {
    fprintf(stderr, "N must be \"Kisuu\". Try again.\n");
    exit(0);
  }
  k = (n + 1) / 2;
  m = n - k;
  
  matrixInit();

  poolptr[1] = pool[1];
  printResult(search(1));
}

void matrixInit()
{
  int i, j;
  int w = 0;
  uShort *pptr;
  
  printf("(n,k,m,d) = (%d,%d,%d,%d)\n", n, k, m, d);
#ifdef NAB
  printf("Select row NAB order ");
#else
  printf("Select row NDB order ");
#endif
#ifdef NAB
  /* column NAB order */
  printf("and column NAB order.\n");
  for (i = exp2(m) - 1, pptr = pool[1], poolnum[1] = 0; i >= 0; i --) {
#else
  /* column NDB order */
  printf("and column NDB order.\n");
  for (i = 0, pptr = pool[1], poolnum[1] = 0; i <= exp2(m) - 1; i ++) {
#endif
    *pptr = 0;
    for (j = 0; j < m; j ++)  {
      weight(i) += ((i & mask(j)) != 0);
    }
    if (weight(i) >= d - 1)  {
      *pptr = i;
      pptr ++;
      poolnum[1] ++;
    }
  }
}

int search(lev)
  int lev;
{
  int i, count = 0;
  uShort u;
  uShort *next;
  int v;
  int orderingFlag = TRUE;
  int number;
  
  do  {
    u = *poolptr[lev];
    for (i = 0; i < m; i ++)  {
      rowValue[lev][i] = rowValue[lev - 1][i] << 1;
      rowValue[lev][i] += ((u & mask(m - i - 1)) != 0);
      if (i
       && orderingFlag
#ifdef NAB
    /* row NAB order */
       && rowValue[lev][i - 1] < rowValue[lev][i])  {
#else
    /* row NDB order */
       && rowValue[lev][i - 1] > rowValue[lev][i])  {
#endif
        orderingFlag = FALSE;
      }
    }
    if (orderingFlag)  {
      if (lev == k)  {
        printResult(TRUE);
        exit(0);
      /* return TRUE;*/
      }
      else  { /* lev < k */
        if (lev < CRITERIA)  {
          number = poolnum[lev] - count;
        }
        else  {
          for (next = poolptr[lev] + 1, number = 0;
               next < &pool[lev][poolnum[lev]];
               next ++)  {
            if (dependTab[lev - 1][*next ^ u] + 1 >= d - 1)  {
              number ++;
            }
          }
        }
        if (number >= k - lev)  {
          poolnum[lev + 1] = number;
          for (v = 0; v < exp2(m); v ++)  {
            dependTab[lev][v] = min(dependTab[lev - 1][v],
              dependTab[lev - 1][v ^ u] + 1);
          }
          for (next = poolptr[lev] + 1, poolnum[lev + 1] = 0;
               next < &pool[lev][poolnum[lev]];
               next ++)  {
            if (dependTab[lev][*next] >= d - 1)  {
              pool[lev + 1][poolnum[lev + 1] ++] = *next;
            }
          }
          poolptr[lev + 1] = &pool[lev + 1][0];
          search(lev + 1);
        }
      }
    }
    orderingFlag = TRUE;
    poolptr[lev] ++;
    count ++;
  } while (poolptr[lev] <= &pool[lev][poolnum[lev] - k + lev]);
  return FALSE;
}

void  printResult(flag)
  int flag;
{
  int i;
  
  if (flag == FALSE)  {
    printf("***No code was found.\n");
  }
  else  {
    printf("***Code was found.  Trans. of Parity check matrix is:\n");
    for (i = m - 1; i >= 0; i --)  {
      printVec(exp2(i));
    }
    for (i = 1; i <= k; i ++)  {
      printVec(*poolptr[i]);
    }
    calWeight();
  }
  printf("\n");
  vtimes(&vtms, NULL);
  printf("***CPU Time = %0.1f + %0.1f.***\n",
    (double)vtms.vm_utime/60,
    (double)vtms.vm_stime/60);
}

void  printVec(v)
  int v;
{
  int i;

  for (i = m - 1; i >= 0; i--)  {
    if (v & exp2(i))  {
      printf("1");
    }
    else  {
      printf("0");
    }
  }
  printf("\n");
}

void  calWeight()
{
  unsigned long weightVec[MAXK], weightOfThisCodeWord;
  int weightCount[MAXK * 2];
  uShort info;
  unsigned int code;
  int i, j;
  
  for (i = 0; i <= n; i ++)  {
    weightCount[i] = 0;
  }

  printf("\n***Code Genetate Matrix is:\n");
  for (i = 0; i < k; i ++)  {
    weightVec[i] = (*poolptr[i + 1] << k) + exp2(k - i - 1);
    for (j = n - 1; j >= 0; j--)  {
      if (weightVec[i] & exp2(j))  {
        printf("1");
      }
      else  {
        printf("0");
      }
    }
    printf("\n");
  }

  printf("***Weights of this code is: (from 0 to %d)\n", n);
  for (info = 0; info < exp2(k); info ++)  {
    code = weightOfThisCodeWord = 0;
    for (i = 0; i < k; i ++)  {
      if ((info & mask(k - i - 1)) != 0)  {
        code ^= weightVec[i];
      }
    }
    for (i = 0; i < n; i ++)  {
      weightOfThisCodeWord += ((code & mask(i)) != 0);
    }
    weightCount[weightOfThisCodeWord] ++;
  }
  
  for (i = 0; i <= n; i ++)  {
    printf("%d ", weightCount[i]);
  }
}
