#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include <limits.h>

#define min(a, b) (a > b ? b : a)
#define max(a, b) (a < b ? b : a)
#define degToRad(angleInDegrees) ((angleInDegrees) * M_PI / 180.0)

const int MAP_SIZE = 21000000;
const int WIDTH = 7000;
const int HEIGHT = 3000;

const int POPULATION_SIZE = 40;
const int CHROMOSOME_SIZE = 10;
const int ACTION_DURATION = 40;
const int ELITISM = 15;

double get_rand() {
  return rand() / (double) RAND_MAX;
}

int get_random(int min, int max) {
  return (int) roundf(get_rand() * (max - min + 1)) + min;;
}

int trim(int value) {
  return value < 2 ? 0 : get_random(1, value);
}

bool random_bool() {
  return get_rand() > 0.5;
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

double distance(Coordinate from, Coordinate to) {
  return sqrt((to.x - from.x) * (to.x - from.x) + (to.y - from.y) * (to.y - from.y));
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
  int duration;
} Action;

typedef struct {
  Action * actions;
  int size;
  double fitness;
} ActionSequence;

int compare_sequences(const void * _a, const void * _b) {
  ActionSequence * a = (ActionSequence *) _a;
  ActionSequence * b = (ActionSequence *) _b;
  if (b->fitness > a->fitness) return -1;
  if (b->fitness < a->fitness) return 1;
  return 0;
}

Action zero_action() {
  return (Action) {0, 0, 0};
}

Action create_random_action() {
  int duration = get_random(1, ACTION_DURATION);
  int rotation = get_random(-90, 90);
  int thrust = get_random(0, 4);
  return (Action) {thrust, rotation, duration};
}

typedef struct {
  ActionSequence * sequences;
  int size;
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

Lander advance(Lander lander, Action action) {
  float xAcc = ((float) lander.action.thrust) * sinf(degToRad((float) lander.action.angle));
  float yAcc = ((float) lander.action.thrust) * cosf(degToRad((float) lander.action.angle)) - 3.711f;

  int nextX = (int) roundf((float) lander.coordinate.x + (float) lander.speed.x - (0.5f * xAcc));
  int nextY = (int) roundf((float) lander.coordinate.y + (float) lander.speed.y + (0.5f * yAcc));

  Coordinate next_coordinate = (Coordinate) {nextX, nextY};
  Vector2D next_speed = (Vector2D) {lander.speed.x - (int) roundf(xAcc), lander.speed.y + (int) roundf(yAcc)};
  Action next_action = (Action) { lander.action.thrust + action.thrust, lander.action.angle + action.angle };
  next_action.thrust = max(0, min(4, next_action.thrust));
  next_action.angle = max(-90, min(90, next_action.angle));

  Lander result = (Lander) {
    next_coordinate,
    next_speed,
    lander.fuel - lander.action.thrust,
    next_action
  };

  return result;
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

bool is_crash(Coordinate coordinate, Surface * surface) {
  int index = coordinate.y * WIDTH + coordinate.x;
  return surface->ground[index];
}

bool is_lost(Coordinate coordinate) {
  return coordinate.x < 0 || coordinate.x >= WIDTH || coordinate.y < 0 || coordinate.y >= HEIGHT;
}

bool is_game_over(Coordinate coordinate, Surface * surface) {
  return is_lost(coordinate) || is_crash(coordinate, surface);
}

bool check_landing_zone(Coordinate coordinate, Surface * surface) {
  int min_x = min(surface->landing[0].x, surface->landing[1].x);
  int max_x = max(surface->landing[0].x, surface->landing[1].x);

  return abs(coordinate.y - surface->landing[0].y) <= 200 && coordinate.x >= min_x && coordinate.x <= max_x;
}

double evaluate_sequence(ActionSequence sequence, Lander lander, Surface * surface) {
  Lander evaluated_lander = lander;
  int current_action = 0;
  int action_duration = 0;

  while (true) {
    if (action_duration == sequence.actions[current_action].duration) {
      current_action += 1;
      action_duration = 0;
    }
    evaluated_lander = advance(evaluated_lander, sequence.actions[current_action]);

    if (is_game_over(evaluated_lander.coordinate, surface)) {
      if (!check_landing_zone(evaluated_lander.coordinate, surface)) return -1.0;
      return 100000.0 + evaluated_lander.fuel -
        (abs(evaluated_lander.speed.y) - 35) * 200 -
        (abs(evaluated_lander.speed.x) - 15) * 200 -
        abs(evaluated_lander.action.angle) * 1000;
    }
    action_duration += 1;
  }
}

void cross_over(ActionSequence a, ActionSequence b, ActionSequence * save_in) {
  ActionSequence start;
  ActionSequence end;

  if (random_bool()) {
    start = b;
    end = a;
  } else {
    start = a;
    end = b;
  }

  int end_start = trim(end.size);
  int start_end = trim(min(start.size, end_start + CHROMOSOME_SIZE - start.size - 1));

  int result_size = start_end + (end.size - end_start);
  save_in->size = result_size;

  for (int i = 0; i < start_end; i++) {
    save_in->actions[i] = start.actions[i];
  }

  for (int i = 0; i < (end.size - end_start); i++) {
    save_in->actions[start_end + i] = end.actions[end_start + i];
  }
}

ActionSequence mutate_sequence(ActionSequence target) {
  return target;
}

Population * initialize_population() {
  // We allocate all the required memory for a population. The initialization results in an empty Population
  Population * population = (Population *) malloc(sizeof(Population));
  population->sequences = (ActionSequence *) malloc(POPULATION_SIZE * sizeof(ActionSequence));
  population->size = 0;
  for (int i = 0; i < POPULATION_SIZE; i++){
    population->sequences[i].actions = (Action *) calloc(CHROMOSOME_SIZE, sizeof(Action));
    population->sequences[i].size = 0;
    population->sequences[i].fitness = -1.0;
  }

  return population;
}

void to_random_sequence(ActionSequence * action_sequence) {
  int size = get_random(1, CHROMOSOME_SIZE);
  action_sequence->size = size;
  action_sequence->fitness = -1.0;
  for (int i = 1; i <= size; i++) {
    if (i == size) {
      action_sequence->actions[i - 1] = (Action) {get_random(0, 4), 0, INT_MAX};
    } else {
      action_sequence->actions[i - 1] = create_random_action();
    }
  }
}

Action move_to_target(Lander lander, Surface * surface) {
  int target_x = (surface->landing[0].x + surface->landing[0].x) / 2;

  bool is_left = lander.coordinate.x > target_x;

  int desired_rotation = is_left ? 10 : -10;

  int rotate = desired_rotation - lander.action.angle;

  return (Action) {1, rotate};
}

int main() {
  srand(time(NULL));
  Surface * surface = create_surface();
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

    double limit = get_millis() + 90;
    Population * population = initialize_population();
    for (int i = 0; i < ELITISM; i++) {
       population->sequences[i].actions[0] = lander.action;
       population->sequences[i].size = 1;
       population->size += 1;
    }
    ActionSequence best;
    best.actions = (Action *) calloc(CHROMOSOME_SIZE, sizeof(Action));
    best.fitness = INT_MIN;
    best.size = 0;

    int iterations = 0;
    while (get_millis() < limit) {
      iterations += 1;
      int population_size = population->size;
      for (int i = population_size; i < POPULATION_SIZE; i++) {
        to_random_sequence(&population->sequences[i]);
        population->size += 1;
      }
      for (int i = 0; i < POPULATION_SIZE; i++) {
        population->sequences[i].fitness = evaluate_sequence(population->sequences[i], lander, surface);
      }
      qsort(population->sequences, POPULATION_SIZE, sizeof(ActionSequence), compare_sequences);
      if (best.fitness < population->sequences[0].fitness) {
        best = population->sequences[0];
      }
      for (int i = 0; i < ELITISM; i++) {
        int change_at = ELITISM + i;
        cross_over(population->sequences[get_random(0, ELITISM)], population->sequences[get_random(0, ELITISM)], &population->sequences[change_at]);
      }
      population->size -= 10; // We will "remove" the last 10 elements, in order to add 10 random actions in the next iteration.
      for (int i = 0; i < population->size; i++) {
        population->sequences[i] = mutate_sequence(population->sequences[i]);
      }
    }

    fprintf(stderr, "Done %d iterations. The best has %f\n", iterations, population->sequences[0].fitness);

    Action better = population->sequences[0].actions[0];

//    if (population->sequences[0].fitness < 0) {
//      better = move_to_target(lander, surface);
//    }

    better.thrust += lander.action.thrust;
    better.angle += lander.action.angle;
    better.thrust = min(4, max(0, better.thrust));
    better.angle = min(90, max(-90, better.angle));

    Lander will_be = advance(lander, better);

    fprintf(stderr, "After action the lander will be in (%d,%d) with s (%d,%d), angle %d and thrust %d\n", will_be.coordinate.x, will_be.coordinate.y, will_be.speed.x, will_be.speed.y, will_be.action.angle, will_be.action.thrust);

    printf("%d %d\n", better.angle, better.thrust);
  }

  return 0;
}
