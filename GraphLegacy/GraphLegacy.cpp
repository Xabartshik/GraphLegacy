// GraphLegacy.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
//Сделал -- Ошлаков Данил, ИВТ-22
#include <iostream>
#include "GraphLegacy.h"

using namespace std;

int main() {
    WeightedGraphLegacy<std::string>::runTests();
    WeightedGraphLegacy<std::string> graph;
    WeightedGraphLegacy<std::string> randGraph;
    randGraph = createRandomGraph(20, 10);
    std::cout << "Сгенерированный: " << std::endl;
    randGraph.printVertices();
    graph.printEdges();
    std::cout << "Количество ребер: " << randGraph.getEdgeCount() << std::endl;
    std::cout << "Количество вершин: " << randGraph.getVertexCount() << std::endl;
    std::cout << "Минимально возможный путь между вершинами А и С: " << std::endl;

    std::cout << "Ручной граф: " << std::endl;

    graph.insertVertex("A");
    graph.insertVertex("B");
    graph.insertVertex("C");
    graph.insertVertex("D");
    graph.insertVertex("E");
    graph.insertEdge("E", "A", 4.0);
    graph.insertEdge("B", "E", 4.0);
    graph.insertEdge("A", "B", 2.0);
    graph.insertEdge("A", "A", 1.0);
    graph.insertEdge("B", "C", 3.0);
    graph.insertEdge("C", "D", 1.0);
    graph.insertEdge("D", "A", 4.0);

    graph.printVertices();
    graph.printEdges();
    std::cout << "Количество ребер: " << graph.getEdgeCount() << std::endl;
    std::cout << "Количество вершин: " << graph.getVertexCount() << std::endl;
    std::cout << "Минимально возможный путь между вершинами А и С: " << std::endl;
    vector<double> distances = graph.findShortestPath("A", "C");
    for (double val : distances)
    {
        cout << val << "   ";
    }
    cout << endl;
    cout << distances.size() << endl;
    graph.depthFirstSearch("A");
    graph.breadthFirstSearch("A");
    graph.printAdjacencyMatrix();
    return 0;
}