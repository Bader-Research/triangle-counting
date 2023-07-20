#include "types.h"
#include "queue.h"

#define EMPTY ((UINT_t) (-1))

// Function to create a new queue
Queue *createQueue(UINT_t size) {
  Queue *queue = (Queue *)malloc(sizeof(Queue));
  assert_malloc(queue);
  queue->items = (UINT_t *)malloc(size * sizeof(UINT_t));
  assert_malloc(queue->items);
  queue->front = EMPTY;
  queue->rear = EMPTY;
  queue->size = size;
  return queue;
}

void free_queue(Queue *queue) {
  free(queue->items);
  free(queue);
}

// Function to check if the queue is empty
int isEmpty(Queue *queue) {
  return queue->rear == EMPTY;
}

// Function to check if the queue is full
int isFull(Queue *queue) {
  return queue->rear == queue->size - 1;
}

// Function to add an element to the queue
void enqueue(Queue *queue, UINT_t value) {
  if (isFull(queue))
    fprintf(stderr,"Queue is full.\n");
  else {
    if (queue->front == EMPTY) {
      queue->front = 0;
      queue->rear  = 0;
      queue->items[queue->rear] = value;
    } else {
      queue->rear++;
      queue->items[queue->rear] = value;
    }
  }
}

// Function to remove an element from the queue
UINT_t dequeue(Queue *queue) {
  UINT_t item;
  if (isEmpty(queue)) {
    printf("Queue is empty.\n");
    item = EMPTY;
  } else {
    item = queue->items[queue->front];
    queue->front++;
    if (queue->front > queue->rear)
      queue->front = queue->rear = EMPTY;
  }
  return item;
}



