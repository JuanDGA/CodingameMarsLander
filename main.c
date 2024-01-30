#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

#define min(a, b) (a > b ? b : a)
#define max(a, b) (a < b ? b : a)

const int MAP_SIZE = 21000000;
const int WIDTH = 7000;
const int HEIGHT = 3000;

const int POPULATION_SIZE = 40;
const int CHROMOSOME_SIZE = 10;
const int ELITISM = 4;

double get_rand() {
  return rand() / (double) RAND_MAX;
}

int get_random(int min, int max) {
  return (int) roundf(get_rand() * (max - min + 1)) + min;;
}

double get_millis() {
  return ((double) clock() / (double) CLOCKS_PER_SEC) * 1000.0;
}

typedef struct {
  int x;
  int y;
} Coordinate;

Coordinate zero_coordinate() {
  return (Coordinate) {0, 0};
}

typedef struct {
  int x;
  int y;
} Vector2D;

Vector2D zero_vector() {
  return (Vector2D) {0, 0};
}

bool is_valid_coordinate(int x, int y) {
  return x >= 0 && x < WIDTH && y >= 0 && y < HEIGHT;
}

typedef struct {
  bool * ground; // 7000 * 3000
  Coordinate landing[2];
} Surface;

Surface * create_surface() {
  Surface * surface = (Surface *) malloc(sizeof(Surface));
  surface->ground = (bool *) calloc(sizeof(bool), MAP_SIZE);
  surface->landing[0].x = 0;
  surface->landing[0].y = 0;
  surface->landing[1].x = 0;
  surface->landing[1].y = 0;
  return surface;
}

typedef struct {
  int thrust;
  int angle;
} Action;

Action zero_action() {
  return (Action) {0, 0};
}

typedef struct {
  Action actions[10];
} Chromosome;

typedef struct {
  Chromosome chromosomes[40];
} Population;

typedef struct {
  Coordinate coordinate;
  Vector2D speed;
  int fuel;
  Action action;
} Lander;

Lander initialize_lander() {
  return (Lander) {zero_coordinate(), zero_vector(), 0, zero_action()};
}

Coordinate advance(Lander lander) {
  float xAcc = ((float) lander.action.thrust) * sinf((float) lander.action.angle);
  float yAcc = ((float) lander.action.thrust) * cosf((float) lander.action.angle) - 3.711f;

  int nextX = (int) roundf((float) lander.coordinate.x + (float) lander.speed.x + (0.5f * xAcc));
  int nextY = (int) roundf((float) lander.coordinate.y + (float) lander.speed.y + (0.5f * yAcc));
  return (Coordinate) {nextX, nextY};
}

void push(int * stack, int * top, int x, int y) {
  stack[++(*top)] = x;
  stack[++(*top)] = y;
}

bool pop(const int * stack, int * top, int * x, int * y) {
  if (*top < 0) return false;
  *y = stack[(*top)--];
  *x = stack[(*top)--];
  return true;
}

void draw_line(Surface * surface, Coordinate from, Coordinate to) {
  int dx = to.x - from.x;
  int dy = to.y - from.y;
  int steps = abs(dx) > abs(dy) ? abs(dx) : abs(dy);

  float xIncrement = (float) dx / (float) steps;
  float yIncrement = (float) dy / (float) steps;

  float x = (float) from.x;
  float y = (float) from.y;

  for (int i = 0; i <= steps; i++) {
    int index = (int) (roundf(y) * (float) WIDTH + roundf(x));
    surface->ground[index] = true;
    x += xIncrement;
    y += yIncrement;
  }
}

void add_point(Surface * surface, Coordinate * coordinates, int current_position) {
  Coordinate next_point = coordinates[current_position];
  int index = next_point.y * WIDTH + next_point.x;
  surface->ground[index] = true;
  if (current_position > 0) {
    Coordinate last_point = coordinates[current_position - 1];
    if (last_point.y == next_point.y) {
      surface->landing[0].x = last_point.x;
      surface->landing[0].y = last_point.y;
      surface->landing[1].x = next_point.x;
      surface->landing[1].y = next_point.y;
    }
    draw_line(surface, last_point, next_point);
  }
}

void flood_fill_scanline(Surface * surface, Coordinate from) {
  int x1;

  int x = from.x;
  int y = from.y;

  int spanAbove, spanBelow;
  int * stack = malloc(42000000 * sizeof(int));

  int top = -1;

  push(stack, &top, x, y);

  while (pop(stack, &top, &x, &y)) {
    x1 = x;
    while (x1 >= 0 && is_valid_coordinate(x1, y) && !surface->ground[y * WIDTH + x1]) x1--;
    x1 += 1;
    spanAbove = spanBelow = 0;

    while (x1 < WIDTH && !surface->ground[y * WIDTH + x1]) {
      surface->ground[y * WIDTH + x1] = true;
      if(!spanAbove && y > 0 && !surface->ground[(y - 1) * WIDTH + x1]) {
        push(stack, &top, x1, y - 1);
        spanAbove = 1;
      }
      else if(spanAbove && y > 0 && surface->ground[(y - 1) * WIDTH + x1]) spanAbove = 0;
      if(!spanBelow && y < HEIGHT - 1 && !surface->ground[(y + 1) * WIDTH + x1]) {
        push(stack, &top, x1, y + 1);
        spanBelow = 1;
      }
      else if(spanBelow && y < HEIGHT - 1 && surface->ground[(y + 1) * WIDTH + x1]) spanBelow = 0;
      x1 += 1;
    }
  }

  free(stack);
}

void fill_surface(Surface * surface) {
  flood_fill_scanline(surface, (Coordinate) { (surface->landing[0].x + surface->landing[1].x) / 2, surface->landing[1].y - 5});
}

Population initialize_population(Population population, bool is_empty) {
  Population new_population;

  int current_size = 0;
  if (!is_empty) {
    new_population = population;
    current_size = ELITISM;
  }

  for (int i = current_size; i < POPULATION_SIZE; i++) {
    Chromosome chromosome = new_population.chromosomes[i];
    for (int j = 0; j < CHROMOSOME_SIZE; j++) {
      chromosome.actions[j] = (Action) { get_random(-1, 1), get_random(-15, 15) };
    }
    new_population.chromosomes[i] = chromosome;
    current_size += 1;
  }

  return new_population;
}

int main() {
  srand(time(NULL));
  Surface *surface = create_surface();
  int surface_n;
  scanf("%d", &surface_n);

  Coordinate coordinates[surface_n + 2];
  coordinates[0] = zero_coordinate();
  add_point(surface, coordinates, 0);
  for (int i = 1; i <= surface_n; i++) {
    int land_x, land_y;
    scanf("%d%d", &land_x, &land_y);
    Coordinate coordinate = (Coordinate) {land_x, land_y};
    coordinates[i] = coordinate;
    add_point(surface, coordinates, i);
  }
  coordinates[surface_n + 1] = (Coordinate) {6999, 0};
  add_point(surface, coordinates, surface_n + 1);
  draw_line(surface, zero_coordinate(), coordinates[surface_n + 1]);

  clock_t start_time, end_time;
  start_time = clock();
  fill_surface(surface);
  end_time = clock();

  fprintf(stderr, "Filled surface in %f ms\n", ((double) (end_time - start_time) / CLOCKS_PER_SEC) * 1000.0);


  // game loop
  while (true) {
    Lander lander = initialize_lander();
    scanf(
        "%d%d%d%d%d%d%d",
        &lander.coordinate.x,
        &lander.coordinate.y,
        &lander.speed.x,
        &lander.speed.y,
        &lander.fuel,
        &lander.action.angle,
        &lander.action.thrust
    );

    double limit = get_millis() + 90; // 90ms

    Population population;
    bool population_is_empty = true;

    int iterations = 0;
    while (get_millis() < limit) {
      population = initialize_population(population, population_is_empty);
      population_is_empty = false;
      iterations += 1;
    }

    fprintf(stderr, "Done %d iterations.\n", iterations);

    Coordinate willBeAt = advance(lander);

    int index = willBeAt.y * WIDTH + willBeAt.x;

    fprintf(stderr, "Will be at (%d, %d) [%d]\n", willBeAt.x, willBeAt.y, index);

    if (willBeAt.x >= WIDTH || willBeAt.x < 0 || willBeAt.y >= HEIGHT || willBeAt.y < 0) {
      fprintf(stderr, "Will be lost!\n");
    } else if (surface->ground[index]) {
      fprintf(stderr, "Will crash!\n");
    }

    Action nextAction = population.chromosomes[0].actions[0];

    nextAction.thrust += lander.action.thrust;
    nextAction.angle += lander.action.angle;
    nextAction.thrust = min(4, max(0, nextAction.thrust));
    nextAction.thrust = min(90, max(-90, nextAction.thrust));

    printf("%d %d\n", nextAction.angle, nextAction.thrust);
  }

  return 0;
}
