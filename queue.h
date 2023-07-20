#ifndef _QUEUE_H
#define _QUEUE_H

typedef struct {
  UINT_t *items;
  UINT_t front;
  UINT_t rear;
  UINT_t size;
} Queue;

Queue *createQueue(UINT_t);
void free_queue(Queue *);
int isEmpty(Queue *);
int isFull(Queue *);
void enqueue(Queue *, UINT_t);
UINT_t dequeue(Queue *);

#endif
