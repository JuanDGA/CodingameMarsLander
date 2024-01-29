#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

#define min(a, b) ((a < b) ? a : b)
#define max(a, b) ((a > b) ? a : b)
#define toRadians(a) (((float) a) * (M_PI / 180.0))
#define distace(a1, a2, b1, b2) sqrtf((a1 * a1) + (a2 * a2))
float position(float initialPosition, float initialSpeed, float acceleration) {
    return initialPosition + initialSpeed + (0.5 * acceleration);
}
float finalSpeed(float initialSpeed, float acceleration, float distance) {
    float value = sqrtf((initialSpeed * initialSpeed) + (2 * acceleration * distance));
    if (distance > 0) return value;
    if (distance < 0) return -value;
    return value;
}

const int SURFACE_WIDTH = 7000;
const int SURFACE_HEIGHT = 3000;
const int CLOCKS_PER_MS = CLOCKS_PER_SEC / 1000;

const int POPULATION_SIZE = 60;
const int CHROMOSOME_SIZE = 15;
const int SELECTION_SIZE = 10;
const int ELITISM_SIZE = 2;
const float MUTATION_RATE = 0.1;

const float gravity = -3.711;

double sinCos[362];

// Surface related structs and functions

typedef struct {
    unsigned short x;
    unsigned short y;
} Point;

typedef struct {
    unsigned char * map;
    int width;
    int height;
    Point lastPoint;
    Point landingPointA;
    Point landingPointB;
    Point landingPoint;
} MarsSurface;

MarsSurface * createSurface(int width, int height) {
    MarsSurface * marsSurface = (MarsSurface * ) calloc(1, sizeof(MarsSurface));

    marsSurface->map = calloc(SURFACE_WIDTH * SURFACE_HEIGHT, sizeof(char));
    marsSurface->width = width;
    marsSurface->height = height;

    return marsSurface;
}

void fillLine(MarsSurface * surface, int x, const int * vector, Point * lastPointPerX, unsigned char * pointsPerX) {
    for (int pointX = min(surface->lastPoint.x, x); pointX < max(surface->lastPoint.x, x); pointX++) {
        float exactPointY = (float)surface->lastPoint.y + (((float)pointX - (float)surface->lastPoint.x) / (float)vector[0]) * (float)vector[1];
        int pointY = (int) roundf(exactPointY);

        surface->map[pointX + pointY * surface->width] = 1;

        Point lastPoint = lastPointPerX[pointX];

        unsigned char amount = pointsPerX[pointX];

        lastPointPerX[pointX] = (Point) {pointX, pointY};

        if (amount % 2 == 0) {
            for (int j = min(pointY, lastPoint.y); j < max(pointY, lastPoint.y); j++) {
                surface->map[pointX + j * surface->width] = 1;
            }
        } else if (pointY < lastPoint.y && surface->map[pointX + pointY * surface->width] == 1) {
            for (int j = 0; j < pointY - 1; j++) {
                surface->map[pointX + j * surface->width] = 0;
            }
            lastPointPerX[pointX] = (Point) {pointX, 0};
        }

        pointsPerX[pointX]++;
    }
}

void addPoint(MarsSurface * surface, int x, int y, Point * lastPointPerX, unsigned char * pointsPerX) {
    int * vector = (int *) calloc(2, sizeof(short));
    vector[0] = x - surface->lastPoint.x;
    vector[1] = y - surface->lastPoint.y;

    if (vector[1] == 0 && abs(vector[0]) >= (SURFACE_WIDTH / 7)) {
        surface->landingPointA = surface->lastPoint;
        surface->landingPointB = (Point) {x, y};
        surface->landingPoint = (Point) {(x + surface->lastPoint.x) / 2, y};
    }

    if (vector[0] != 0) fillLine(surface, x, vector, lastPointPerX, pointsPerX);

    surface->lastPoint.x = x;
    surface->lastPoint.y = y;

    free(vector);
}

void printSurface(MarsSurface * surface) {
    printf("Printing...\n");
    for (int y = 0; y < SURFACE_HEIGHT; y++) {
        for (int x = 0; x < SURFACE_WIDTH; x++) {
            int index = x + ((SURFACE_HEIGHT - y - 1) * SURFACE_WIDTH);
            printf("%d", surface->map[index]);
        }
        printf("\n");
    }
}

// Genetic algorithm related structs and functions

typedef struct {
    float xPos, yPos, xSpeed, ySpeed, fuel, angle, thrust;
} Ship;

typedef struct {
    int rotation;
    int thrustChange;
} Action;

typedef struct {
    Action * actions;
    float fitness;
} Chromosome;

typedef struct {
    Chromosome * chromosomes;
} Population;

void printShip(Ship * ship) {
    fprintf(stderr, "Ship: {X=%.2f Y=%.2f HSpeed=%.2f VSpeed=%.2f Fuel=%.2f Angle=%.2f Thrust=%.2f}\n",
            ship->xPos, ship->yPos, ship->xSpeed, ship->ySpeed, ship->fuel, ship->angle, ship->thrust);
}

Ship advanceShip(const Ship * ship, Action action) {
    Ship nextShip = * ship;

    float xAcceleration = nextShip.thrust * sinCos[90 + ((int) nextShip.angle)];
    float yAcceleration = nextShip.thrust * sinCos[182 + 90 + ((int) nextShip.angle)] + gravity;

    float newXPosition = position(nextShip.xPos, nextShip.xSpeed, xAcceleration);
    float newYPosition = position(nextShip.yPos, nextShip.ySpeed, yAcceleration);

    float newXSpeed = finalSpeed(nextShip.xSpeed, xAcceleration, newXPosition - nextShip.xPos);
    float newYSpeed = finalSpeed(nextShip.ySpeed, yAcceleration, newYPosition - nextShip.yPos);

    nextShip.xPos = roundf(newXPosition);
    nextShip.yPos = roundf(newYPosition);
    nextShip.xSpeed = roundf(newXSpeed);
    nextShip.ySpeed = roundf(newYSpeed);
    nextShip.fuel -= nextShip.thrust;
    nextShip.angle += action.rotation;
    nextShip.thrust += action.thrustChange;
    nextShip.angle = fmin(90.0f, fmax(-90.0f, nextShip.angle));
    nextShip.thrust = fmin(4.0f, fmax(0.0f, nextShip.thrust));

    return nextShip;
}

float distance(Point p1, Point p2) {
    float dx = (float)(p2.x - p1.x);
    float dy = (float)(p2.y - p1.y);
    return sqrtf((dx * dx) + (dy * dy));
}

float getAngle(float a, float b, float c) {
    float cosA = ((b * b) + (c * c) - (a * a)) / (2.0 * b * c);
    float angle = acosf(cosA);
    angle = angle * 180.0 / M_PI;
    return angle;
}

float calculateFitness(MarsSurface * surface, Ship * ship) {
    Point shipPoint = (Point) {(int) ship->xPos, (int) ship->yPos};
    float distanceToBase = distance(shipPoint, surface->landingPoint);
    float b = distance(shipPoint, (Point) {surface->landingPoint.x, 0});
    float c = surface->landingPoint.y;

    int minLandingX = min(surface->landingPointA.x, surface->landingPointB.x);
    int maxLandingX = max(surface->landingPointA.x, surface->landingPointB.x);

    float relativeAngle = fabsf(getAngle(distanceToBase, b, c));

    int isLost = ship->xPos >= SURFACE_WIDTH || ship->xPos < 0 || ship->yPos >= SURFACE_HEIGHT;

    int position;

    if (isLost) {
        position = 0;
    } else {
        position = surface->map[((int) ship->xPos) + ((int) ship->yPos) * surface->width];
    }

    int isLandingZone = ship->xPos >= minLandingX &&
            ship->xPos <= maxLandingX &&
            ship->yPos >= surface->landingPoint.y - 20 &&
            ship->yPos <= surface->landingPoint.y;

    int isCrash = ((int) ship->angle) != 0 || abs((int) ship->xSpeed) > 20 || abs((int) ship->ySpeed) > 40;

    int isLandingCrash = isLandingZone && isCrash && position == 1;

    int isCommonCrash = isLandingZone != 1 && position == 1;

    int isFlying = position == 0;

    int isLanding = isLandingZone && isCrash == 0 && position == 1;

    // TODO: Implement Better Fitness calculation

    // CASES:
    // Landing crash -> 80K - speed - angle
    // Common crash -> 30k - distance - relative_angle * intersections
    // Lost -> 30K - distance - relative_angle * intersections
    // Flying -> 50k - distance - relative_angle * intersections - speed - angle
    // Landing -> 100K + Fuel

    float fitness = 0.0;

    if (isLandingCrash) {
        fitness = 80000 - (abs((int) ship->xSpeed) * 200) - (abs((int) ship->ySpeed) * 100) - (abs((int) ship->angle) * 1000);
    } else if (isCommonCrash || isLost) {
        fitness = - distanceToBase * 10 - relativeAngle * 5;
    } else if (isFlying) {
        fitness = 50000 - distanceToBase - relativeAngle * 5 - (abs((int) ship->xSpeed) * 200) - (abs((int) ship->ySpeed) * 100);
    } else if (isLanding) {
        fitness = 100000 + ship->fuel * 10;
    }

    return fitness;

//    return 100000.0 + 2 * ship->fuel -
//        max(0.0, abs((int) ship->xSpeed) - 45) * 200 -
//        max(0.0, abs((int) ship->ySpeed) - 25) * 200 -
//        abs((int) ship->angle) * 1000 -
//        a * 5 - fabsf(getAngle(a, b, c)) * 10;
}

int compareChromosomes(const void * a, const void * b) {
    Chromosome * chromosomeA = (Chromosome *) a;
    Chromosome * chromosomeB = (Chromosome *) b;

    return chromosomeB->fitness - chromosomeA->fitness;
}

Action generateRandomAction() {
    int rotation = (rand() % 31) - 15;
    int thrust = (rand() % 3) - 1;

    return (Action) {rotation, thrust};
}

void initializeChromosome(MarsSurface * surface, Chromosome * chromosome, Ship * ship) {
    chromosome->actions = (Action *) malloc(CHROMOSOME_SIZE * sizeof(Action));

    Ship currentShip = * ship;
    int finished = 0;
    for (int j = 0; j < CHROMOSOME_SIZE; j++) {
        Action action = generateRandomAction();
        chromosome->actions[j] = action;

        if (currentShip.xPos < 0 || currentShip.xPos >= SURFACE_WIDTH || currentShip.yPos >= SURFACE_HEIGHT) finished = 1;
        else if (surface->map[((int) currentShip.xPos) + ((int) currentShip.yPos) * SURFACE_WIDTH]) finished = 1;

        if (finished == 0) {
            currentShip = advanceShip(&currentShip, action);
        }
    }
    chromosome->fitness = calculateFitness(surface, &currentShip);
}

void initializePopulation(MarsSurface * surface, Population * population, Ship * ship) {
    for (int i = 0; i < POPULATION_SIZE; i++) {
        initializeChromosome(surface, &population->chromosomes[i], ship);
    }
}

Population * selectPopulation(Population * basePopulation, Population * selection, Population * tournament) {
    for (int i = 0; i < ELITISM_SIZE; i++) {
        selection->chromosomes[i] = basePopulation->chromosomes[i];
    }

    for (int i = ELITISM_SIZE; i < SELECTION_SIZE; i++) {
        for (int j = 0; j < 5; j++) {
            int idx = rand() % POPULATION_SIZE;
            tournament->chromosomes[j] = basePopulation->chromosomes[idx];
        }
        qsort(tournament->chromosomes, 5, sizeof(Chromosome), compareChromosomes);
        selection->chromosomes[i] = tournament->chromosomes[0];
    }
    return selection;
}

void crossoverSelection(MarsSurface * surface, Population * population, Population * selection, Ship * ship) {
    for (int i = 0; i < POPULATION_SIZE; i++) {
        int selectedA = rand() % 2;
        int selectedB = rand() % 2;

        int splitIdx = rand() % CHROMOSOME_SIZE;

        Ship currentShip = *ship;
        int finished = 0;
        for (int j = 0; j < CHROMOSOME_SIZE; j++) {
            Action selectedAction;

            if (j < splitIdx) {
                selectedAction = selection->chromosomes[selectedA].actions[j];
            } else {
                selectedAction = selection->chromosomes[selectedB].actions[j];
            }

            if ((((float) rand()) / RAND_MAX) < MUTATION_RATE) {
                selectedAction.thrustChange += ((rand() % 3) - 1);
                selectedAction.thrustChange = min(1, max(-1, selectedAction.thrustChange));
                selectedAction.rotation += ((rand() % 31) - 15);
                selectedAction.rotation = min(15, max(-15, selectedAction.rotation));
            }

            if (finished == 0) {
                currentShip = advanceShip(&currentShip, selectedAction);
                if (currentShip.xPos < 0 || currentShip.xPos >= SURFACE_WIDTH || currentShip.yPos >= SURFACE_HEIGHT) finished = 1;
                else if (surface->map[((int) currentShip.xPos) + ((int) currentShip.yPos) * SURFACE_WIDTH]) finished = 1;
            }

            population->chromosomes[i].actions[j] = selectedAction;
        }
        population->chromosomes[i].fitness = calculateFitness(surface, &currentShip);
    }
}

void initializeTrignometricValues() {
    for (int i = -90; i <= 90; i++) {
        int sinPos = i + 90;
        int cosPos = i + 90 + 181;

        float inRadians = toRadians(i);

        sinCos[sinPos] = sinf(inRadians);
        sinCos[sinPos] = cosf(inRadians);
    }
}

int main() {
    srand(time(NULL));
    initializeTrignometricValues();

    MarsSurface * marsSurface = createSurface(SURFACE_WIDTH, SURFACE_HEIGHT);

    Point * lastPointPerX = calloc(SURFACE_WIDTH, sizeof(Point));

    unsigned char * pointsPerX = calloc(SURFACE_WIDTH, sizeof(char));

    Ship * state = (Ship *) malloc(sizeof(Ship));
    float xPos, yPos, hSpeed, vSpeed, remainingFuel, angle, thrust;

    Population * population = (Population *) malloc(sizeof(Population));
    population->chromosomes = malloc(POPULATION_SIZE * sizeof(Chromosome));

    Population * selection = (Population *) malloc(sizeof(Population));
    selection->chromosomes = malloc(SELECTION_SIZE * sizeof(Chromosome));

    Population * tournament = (Population *) malloc(sizeof(Population));
    tournament->chromosomes = malloc(5 * sizeof(Chromosome));

    Action action;

    int vertices, x, y;
    scanf("%d", &vertices);
    clock_t startTime = clock();
    for (int i = 0; i < vertices; i++) {
        scanf("%d%d", &x, &y);
//        fprintf(stderr, "addPoint(marsSurface, %d, %d, lastPointPerX, pointsPerX);\n", x / 100, y / 100);
        addPoint(marsSurface, x, y, lastPointPerX, pointsPerX);
    }

    free(lastPointPerX);
    free(pointsPerX);

    fprintf(stderr, "Landing Point: %d_%d\n", marsSurface->landingPoint.x, marsSurface->landingPoint.y);

//    printSurface(marsSurface);
    double maxDuration = 0.0900 * CLOCKS_PER_SEC;
    short turn = 0;
    while (1) {
        turn += 1;

        if (turn != 1) {
            startTime = clock();
            maxDuration = 0.090 * CLOCKS_PER_SEC;
        }

        scanf("%f%f%f%f%f%f%f",
              &state->xPos, &state->yPos, &state->xSpeed, &state->ySpeed, &state->fuel, &state->angle, &state->thrust
        );
        printShip(state);

        fprintf(stderr, "Took: %ldms\n", (clock() - startTime) / 1000);

        if (turn == 1) initializePopulation(marsSurface, population, state);

        fprintf(stderr, "After initialization: %ldms\n", (clock() - startTime) / 1000);

        int iterations = 0;

        while ((clock() - startTime) < maxDuration) {
            qsort(population->chromosomes, POPULATION_SIZE, sizeof(Chromosome), compareChromosomes);
            crossoverSelection(marsSurface, population, selectPopulation(population, selection, tournament), state);
            iterations++;
        }
        qsort(population->chromosomes, POPULATION_SIZE, sizeof(Chromosome), compareChromosomes);

        fprintf(stderr, "Finished iterations after: %ldms\n", (clock() - startTime) / 1000);
        fprintf(stderr, "Generations: %d\n", iterations);

        action = population->chromosomes->actions[0];

        Ship result = advanceShip(state, action);
        printShip(&result);
        float resultFitness = calculateFitness(marsSurface, &result);
        fprintf(stderr, "Fitness: %f\n", resultFitness);

        printf("%d %d\n", (int) result.angle, (int) result.thrust);
    }

    return 1;
}