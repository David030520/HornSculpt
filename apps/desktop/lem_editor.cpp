#include "lem_editor.hpp"
#include "lem_driver_presets.hpp"
#include "numeric_input.hpp"

#include <QApplication>
#include <QCheckBox>
#include <QComboBox>
#include <QDrag>
#include <QDragEnterEvent>
#include <QDragMoveEvent>
#include <QDropEvent>
#include <QFormLayout>
#include <QGraphicsEllipseItem>
#include <QGraphicsPathItem>
#include <QGraphicsScene>
#include <QGraphicsSceneMouseEvent>
#include <QGraphicsTextItem>
#include <QGraphicsView>
#include <QHBoxLayout>
#include <QHeaderView>
#include <QKeyEvent>
#include <QLabel>
#include <QLineEdit>
#include <QListWidget>
#include <QLocale>
#include <QMimeData>
#include <QPainter>
#include <QPushButton>
#include <QScrollArea>
#include <QSignalBlocker>
#include <QSplitter>
#include <QTableWidget>
#include <QTimer>
#include <QVBoxLayout>
#include <QWheelEvent>
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <map>
#include <numbers>
#include <set>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace hornsculpt {
namespace {
constexpr auto componentMime = "application/x-hornsculpt-lem-component";
constexpr int itemKindRole = 0;
constexpr int itemIndexRole = 1;
constexpr int terminalRole = 2;
enum class DiagramItemKind { component = 1, node, wire };
bool portFace(const LemFaceItem &face, lem::BemPortKind kind) {
    return kind == lem::BemPortKind::interfacePort
               ? face.boundary == BoundaryCondition::wall && face.interfaceAllowed
               : face.boundary == BoundaryCondition::driven;
}

void retainPortWeights(lem::BemPort &port) {
    std::erase_if(port.velocityWeights, [&](const auto &entry) {
        return std::find(port.faceIds.begin(), port.faceIds.end(), entry.first) ==
               port.faceIds.end();
    });
}

QString portWeightsText(const std::map<std::uint32_t, double> &weights) {
    QStringList entries;
    for (const auto &[id, weight] : weights)
        entries.push_back(QString::number(id) + ':' + QString::number(weight, 'g', 17));
    return entries.join("; ");
}

std::map<std::uint32_t, double> parsedPortWeights(const QString &input, const lem::BemPort &port) {
    std::map<std::uint32_t, double> weights;
    if (input.trimmed().isEmpty())
        return weights;
    for (const auto &entry : input.split(';')) {
        const auto pair = entry.split(':');
        if (pair.size() != 2)
            throw std::invalid_argument(
                "Format: Flächen-ID:Faktor; weitere Einträge mit Semikolon.");
        bool idValid = false, weightValid = false;
        const auto id = pair[0].trimmed().toUInt(&idValid);
        const auto weight = QLocale::c().toDouble(pair[1].trimmed(), &weightValid);
        if (!idValid || !id ||
            std::find(port.faceIds.begin(), port.faceIds.end(), id) == port.faceIds.end())
            throw std::invalid_argument("Gewichte dürfen nur zugeordnete Flächen-IDs enthalten.");
        if (!weightValid || !std::isfinite(weight) || weight <= 0)
            throw std::invalid_argument(
                "Gewicht muss endlich und größer als 0 sein; Dezimalpunkt verwenden.");
        if (!weights.emplace(id, weight).second)
            throw std::invalid_argument("Jede Flächen-ID darf nur ein Gewicht besitzen.");
    }
    return weights;
}

struct TerminalAddress {
    std::size_t element{};
    std::size_t terminal{};
    bool operator==(const TerminalAddress &) const = default;
};

QColor domainColor(lem::Domain domain) {
    switch (domain) {
    case lem::Domain::electrical:
        return QColor(180, 85, 35);
    case lem::Domain::mechanical:
        return QColor(70, 135, 80);
    case lem::Domain::acoustic:
        return QColor(50, 115, 185);
    }
    return QApplication::palette().color(QPalette::Text);
}

class ComponentPalette final : public QListWidget {
  public:
    explicit ComponentPalette(QWidget *parent = nullptr) : QListWidget(parent) {
        setDragEnabled(true);
        setSelectionMode(QAbstractItemView::SingleSelection);
    }

  protected:
    void startDrag(Qt::DropActions) override {
        if (!currentItem())
            return;
        auto *mime = new QMimeData;
        mime->setData(componentMime, currentItem()->data(Qt::UserRole).toByteArray());
        auto *drag = new QDrag(this);
        drag->setMimeData(mime);
        drag->exec(Qt::CopyAction);
    }
};

class CircuitView final : public QGraphicsView {
  public:
    explicit CircuitView(QGraphicsScene *scene, QWidget *parent = nullptr)
        : QGraphicsView(scene, parent) {
        setAcceptDrops(true);
        setDragMode(QGraphicsView::RubberBandDrag);
        setRenderHint(QPainter::Antialiasing);
        setTransformationAnchor(QGraphicsView::AnchorUnderMouse);
        setMinimumSize(360, 240);
    }
    std::function<void(QByteArray, QPointF)> onDrop;
    std::function<void()> onDelete;
    std::function<void()> onCancelWire;

  protected:
    void dragEnterEvent(QDragEnterEvent *event) override {
        if (event->mimeData()->hasFormat(componentMime))
            event->acceptProposedAction();
        else
            QGraphicsView::dragEnterEvent(event);
    }
    void dragMoveEvent(QDragMoveEvent *event) override {
        if (event->mimeData()->hasFormat(componentMime))
            event->acceptProposedAction();
        else
            QGraphicsView::dragMoveEvent(event);
    }
    void dropEvent(QDropEvent *event) override {
        if (event->mimeData()->hasFormat(componentMime)) {
            if (onDrop)
                onDrop(event->mimeData()->data(componentMime),
                       mapToScene(event->position().toPoint()));
            event->acceptProposedAction();
        } else {
            QGraphicsView::dropEvent(event);
        }
    }
    void keyPressEvent(QKeyEvent *event) override {
        if (event->key() == Qt::Key_Delete && onDelete) {
            onDelete();
            event->accept();
        } else if (event->key() == Qt::Key_Escape && onCancelWire) {
            onCancelWire();
            event->accept();
        } else {
            QGraphicsView::keyPressEvent(event);
        }
    }
    void wheelEvent(QWheelEvent *event) override {
        if (event->modifiers().testFlag(Qt::ControlModifier)) {
            const double factor = event->angleDelta().y() > 0 ? 1.12 : 1 / 1.12;
            const auto current = transform().m11();
            if (current > 0) {
                const auto next = std::clamp(current * factor, .25, 3.0);
                scale(next / current, next / current);
            }
            event->accept();
        } else {
            QGraphicsView::wheelEvent(event);
        }
    }
};

class ComponentItem final : public QGraphicsRectItem {
  public:
    ComponentItem(std::size_t index, const QString &title, const QString &subtitle, double height)
        : QGraphicsRectItem(0, 0, 200, height) {
        setData(itemKindRole, static_cast<int>(DiagramItemKind::component));
        setData(itemIndexRole, QVariant::fromValue<qulonglong>(index));
        setFlags(ItemIsMovable | ItemIsSelectable | ItemSendsGeometryChanges);
        setBrush(QApplication::palette().brush(QPalette::Base));
        setPen(QPen(QApplication::palette().color(QPalette::Mid), 1));
        auto *heading = new QGraphicsTextItem(title, this);
        heading->setDefaultTextColor(QApplication::palette().color(QPalette::Text));
        heading->setPos(12, 5);
        heading->setTextWidth(176);
        auto *detail = new QGraphicsTextItem(subtitle, this);
        detail->setDefaultTextColor(QApplication::palette().color(QPalette::Text));
        detail->setPos(12, 29);
        detail->setTextWidth(176);
    }
    std::function<void()> onMoving;
    std::function<void(QPointF)> onMoved;

  protected:
    QVariant itemChange(GraphicsItemChange change, const QVariant &value) override {
        auto result = QGraphicsRectItem::itemChange(change, value);
        if (change == ItemPositionHasChanged && onMoving)
            onMoving();
        return result;
    }
    void mousePressEvent(QGraphicsSceneMouseEvent *event) override {
        initialPosition_ = pos();
        QGraphicsRectItem::mousePressEvent(event);
    }
    void mouseReleaseEvent(QGraphicsSceneMouseEvent *event) override {
        QGraphicsRectItem::mouseReleaseEvent(event);
        if (pos() != initialPosition_ && onMoved)
            onMoved(pos());
    }

  private:
    QPointF initialPosition_;
};

class TerminalItem final : public QGraphicsEllipseItem {
  public:
    TerminalItem(QPointF position, QColor color, QString label, QGraphicsItem *parent)
        : QGraphicsEllipseItem(-7, -7, 14, 14, parent) {
        setPos(position);
        setBrush(color);
        setPen(QPen(QApplication::palette().color(QPalette::Base), 1));
        setToolTip(std::move(label));
        setCursor(Qt::CrossCursor);
        setZValue(3);
    }
    std::function<void()> onClicked;

  protected:
    void mousePressEvent(QGraphicsSceneMouseEvent *event) override {
        event->accept();
    }
    void mouseReleaseEvent(QGraphicsSceneMouseEvent *event) override {
        event->accept();
        if (onClicked)
            onClicked();
    }
};

class JunctionItem final : public QGraphicsEllipseItem {
  public:
    JunctionItem(std::size_t index, bool reference, const QString &label, QColor color)
        : QGraphicsEllipseItem(-5, -5, 10, 10) {
        setFlags(ItemIsMovable | ItemIsSelectable | ItemSendsGeometryChanges);
        setData(itemKindRole, static_cast<int>(DiagramItemKind::node));
        setData(itemIndexRole, QVariant::fromValue<qulonglong>(index));
        setBrush(color);
        setPen(QPen(color, 1));
        setZValue(2);
        auto *text = new QGraphicsTextItem(label, this);
        text->setDefaultTextColor(QApplication::palette().color(QPalette::Text));
        text->setPos(8, -14);
        if (reference) {
            auto *ground = new QGraphicsPathItem(this);
            QPainterPath path;
            path.moveTo(0, 5);
            path.lineTo(0, 13);
            path.moveTo(-12, 13);
            path.lineTo(12, 13);
            path.moveTo(-8, 18);
            path.lineTo(8, 18);
            path.moveTo(-4, 23);
            path.lineTo(4, 23);
            ground->setPath(path);
            ground->setPen(QPen(color, 1.5));
        }
    }
    std::function<void()> onMoving;
    std::function<void(QPointF)> onMoved;
    std::function<void()> onClicked;

  protected:
    QVariant itemChange(GraphicsItemChange change, const QVariant &value) override {
        auto result = QGraphicsEllipseItem::itemChange(change, value);
        if (change == ItemPositionHasChanged && onMoving)
            onMoving();
        return result;
    }
    void mousePressEvent(QGraphicsSceneMouseEvent *event) override {
        initialPosition_ = pos();
        QGraphicsEllipseItem::mousePressEvent(event);
    }
    void mouseReleaseEvent(QGraphicsSceneMouseEvent *event) override {
        QGraphicsEllipseItem::mouseReleaseEvent(event);
        if (pos() != initialPosition_ && onMoved)
            onMoved(pos());
        else if (onClicked)
            onClicked();
    }

  private:
    QPointF initialPosition_;
};

QString complexText(std::complex<double> value, int precision = 4) {
    return QStringLiteral("%1 %2 j%3")
        .arg(value.real(), 0, 'g', precision)
        .arg(value.imag() < 0 ? QStringLiteral("−") : QStringLiteral("+"))
        .arg(std::abs(value.imag()), 0, 'g', precision);
}

QString domainName(lem::Domain domain) {
    switch (domain) {
    case lem::Domain::electrical:
        return QObject::tr("Elektrisch");
    case lem::Domain::mechanical:
        return QObject::tr("Mechanisch");
    case lem::Domain::acoustic:
        return QObject::tr("Akustisch");
    }
    return {};
}

struct Terminal {
    QString label;
    lem::NodeId *node{};
    std::optional<lem::Domain> domain;
};

std::vector<Terminal> terminals(lem::Element &element) {
    return std::visit(
        [](auto &component) -> std::vector<Terminal> {
            using T = std::decay_t<decltype(component)>;
            if constexpr (std::is_same_v<T, lem::Driver>) {
                return {{QObject::tr("Elektrisch +"), &component.electrical.positive,
                         lem::Domain::electrical},
                        {QObject::tr("Elektrisch −"), &component.electrical.negative,
                         lem::Domain::electrical},
                        {QObject::tr("Vorderseite"), &component.front, lem::Domain::acoustic},
                        {QObject::tr("Rückseite"), &component.rear, lem::Domain::acoustic}};
            } else {
                const auto domain = []() -> std::optional<lem::Domain> {
                    if constexpr (std::is_same_v<T, lem::Cavity> ||
                                  std::is_same_v<T, lem::BemPort> ||
                                  std::is_same_v<T, lem::AcousticSlit> ||
                                  std::is_same_v<T, lem::ClosedSlitCavity>)
                        return lem::Domain::acoustic;
                    return std::nullopt;
                }();
                return {{QObject::tr("Anschluss +"), &component.terminals.positive, domain},
                        {QObject::tr("Anschluss −"), &component.terminals.negative, domain}};
            }
        },
        element.component);
}

QString componentName(const lem::Component &component) {
    return std::visit(
        [](const auto &value) -> QString {
            using T = std::decay_t<decltype(value)>;
            if constexpr (std::is_same_v<T, lem::Resistance>)
                return QObject::tr("Widerstand");
            else if constexpr (std::is_same_v<T, lem::Inertance>)
                return QObject::tr("Induktivität / Masse");
            else if constexpr (std::is_same_v<T, lem::Compliance>)
                return QObject::tr("Kapazität / Nachgiebigkeit");
            else if constexpr (std::is_same_v<T, lem::EffortSource>)
                return QObject::tr("Potentialquelle");
            else if constexpr (std::is_same_v<T, lem::FlowSource>)
                return QObject::tr("Flussquelle");
            else if constexpr (std::is_same_v<T, lem::Cavity>)
                return QObject::tr("Luftkammer");
            else if constexpr (std::is_same_v<T, lem::AcousticSlit>)
                return QObject::tr("Spaltkanal");
            else if constexpr (std::is_same_v<T, lem::ClosedSlitCavity>)
                return QObject::tr("Ringkammer (1D)");
            else if constexpr (std::is_same_v<T, lem::Driver>)
                return QObject::tr("Elektrodynamischer Treiber");
            else if constexpr (std::is_same_v<T, lem::BemPort>)
                return value.kind == lem::BemPortKind::interfacePort
                           ? QObject::tr("BEM↔LEM-Interface")
                           : QObject::tr("BEM-Antriebsanschluss");
        },
        component);
}

struct PaletteEntry {
    const char *id;
    const char *label;
};
constexpr std::array<PaletteEntry, 24> paletteEntries{{{"driver", "Treiber"},
                                                       {"cavity", "Luftkammer"},
                                                       {"a-slit", "Spaltkanal"},
                                                       {"a-closed-slit", "Ringkammer (1D)"},
                                                       {"bem", "BEM-Anschluss"},
                                                       {"bem-interface", "BEM↔LEM-Interface"},
                                                       {"e-source", "Spannungsquelle"},
                                                       {"e-flow", "Stromquelle"},
                                                       {"e-r", "Widerstand [Ω]"},
                                                       {"e-l", "Induktivität [H]"},
                                                       {"e-c", "Kapazität [F]"},
                                                       {"e-ground", "Elektrischer Bezug"},
                                                       {"a-source", "Druckquelle"},
                                                       {"a-flow", "Volumenstromquelle"},
                                                       {"a-r", "Akustischer Widerstand"},
                                                       {"a-l", "Akustische Masse"},
                                                       {"a-c", "Akustische Nachgiebigkeit"},
                                                       {"a-ground", "Druckbezug"},
                                                       {"m-source", "Kraftquelle"},
                                                       {"m-flow", "Geschwindigkeitsquelle"},
                                                       {"m-r", "Mechanischer Widerstand"},
                                                       {"m-l", "Bewegte Masse"},
                                                       {"m-c", "Mechanische Nachgiebigkeit"},
                                                       {"m-ground", "Mechanischer Bezug"}}};
} // namespace

struct LemEditor::Impl {
    explicit Impl(LemEditor &editor) : owner(editor) {}
    LemEditor &owner;
    lem::Network network;
    std::vector<LemFaceItem> faces;
    std::vector<std::uint32_t> selectedFaces;
    std::optional<lem::Solution> solution;
    bool stale{};
    bool rebuilding{};
    std::optional<TerminalAddress> pendingTerminal;
    std::optional<lem::ElementId> selectedElement;
    std::optional<lem::NodeId> selectedNode;
    QCheckBox *enabled{};
    ComponentPalette *palette{};
    QGraphicsScene *scene{};
    CircuitView *view{};
    QScrollArea *properties{};
    QLabel *status{};
    QLabel *resultStatus{};
    QTableWidget *results{};
    QTableWidget *portResults{};
    std::vector<ComponentItem *> componentItems;
    std::vector<JunctionItem *> nodeItems;
    struct Wire {
        TerminalAddress address;
        QGraphicsPathItem *item{};
        TerminalItem *terminal{};
    };
    std::vector<Wire> wires;

    void createUi();
    void rebuildScene();
    void updateWires();
    void commitPositions();
    void updateProperties();
    void updateStatus();
    void updateResults();
    void selectionChanged();
    void notify(bool rebuild = true);
    void addComponent(const QByteArray &kind, QPointF position);
    void addDriverCircuit();
    void addInterfaceCircuit();
    void deleteSelection();
    void terminalClicked(TerminalAddress address);
    void joinNodes(lem::NodeId from, lem::NodeId to);
    void attachTerminal(TerminalAddress address, lem::NodeId target);
    void removeUnusedNodes();
    lem::NodeId createNode(lem::Domain domain, QPointF position, bool reference = false);
    lem::Node *node(lem::NodeId id);
    lem::Element *element(lem::ElementId id);
    lem::Domain terminalDomain(TerminalAddress address);
    lem::NodeId *terminalNode(TerminalAddress address);
    void showError(const QString &text);
};

lem::Node *LemEditor::Impl::node(lem::NodeId id) {
    const auto it = std::find_if(network.nodes.begin(), network.nodes.end(),
                                 [id](const auto &value) { return value.id == id; });
    return it == network.nodes.end() ? nullptr : &*it;
}

lem::Element *LemEditor::Impl::element(lem::ElementId id) {
    const auto it = std::find_if(network.elements.begin(), network.elements.end(),
                                 [id](const auto &value) { return value.id == id; });
    return it == network.elements.end() ? nullptr : &*it;
}

lem::NodeId *LemEditor::Impl::terminalNode(TerminalAddress address) {
    if (address.element >= network.elements.size())
        return nullptr;
    auto pins = terminals(network.elements[address.element]);
    return address.terminal < pins.size() ? pins[address.terminal].node : nullptr;
}

lem::Domain LemEditor::Impl::terminalDomain(TerminalAddress address) {
    auto pins = terminals(network.elements.at(address.element));
    const auto &pin = pins.at(address.terminal);
    if (pin.domain)
        return *pin.domain;
    if (const auto *value = node(*pin.node))
        return value->domain;
    return lem::Domain::electrical;
}

lem::NodeId LemEditor::Impl::createNode(lem::Domain domain, QPointF position, bool reference) {
    std::uint32_t next = 1;
    for (const auto &value : network.nodes)
        next = std::max(next, value.id.value + 1);
    const lem::NodeId id{next};
    network.nodes.push_back(
        {id,
         domain,
         reference,
         (reference ? QStringLiteral("0_%1") : QStringLiteral("n%1")).arg(next).toStdString(),
         {position.x(), position.y()}});
    return id;
}

void LemEditor::Impl::createUi() {
    auto *layout = new QVBoxLayout(&owner);
    layout->setContentsMargins(6, 6, 6, 6);
    auto *toolbar = new QHBoxLayout;
    enabled = new QCheckBox(LemEditor::tr("BEM↔LEM-Kopplung aktiv"));
    enabled->setObjectName(QStringLiteral("lemEnabled"));
    enabled->setToolTip(
        LemEditor::tr("Berechnet das elektrische und akustische Netzwerk mit "
                      "Druckrückwirkung der BEM-Flächen. Ein unvollständiger "
                      "Schaltplan kann bei ausgeschalteter Kopplung bearbeitet werden."));
    enabled->setWhatsThis(enabled->toolTip());
    toolbar->addWidget(enabled);
    auto *preset = new QPushButton(LemEditor::tr("Treiberkreis hinzufügen"));
    preset->setObjectName(QStringLiteral("lemAddDriverCircuit"));
    preset->setToolTip(LemEditor::tr("Fügt einen editierbaren TI-100-Parametersatz, 27 V RMS, "
                                     "0,34 l Rückkammer und einen BEM-Anschluss hinzu. "
                                     "Mmd = 4 g ist eine Annahme aus Hersteller-Mms und muss "
                                     "ohne externe Luftlast verifiziert oder ersetzt werden. "
                                     "Die Spannung ist keine Zusage einer bestimmten Leistung."));
    preset->setWhatsThis(preset->toolTip());
    toolbar->addWidget(preset);
    auto *interface = new QPushButton(LemEditor::tr("Interface + Luftkammer"));
    interface->setObjectName("lemAddInterfaceCircuit");
    interface->setToolTip(
        LemEditor::tr("Verknüpft ausgewählte virtuelle Wall-Flächen mit einem BEM↔LEM-Interface "
                      "und einer Luftkammer. Das Volumen bleibt ungesetzt; Berechnung erst nach "
                      "Eingabe. Keine FEM- oder BEM↔BEM-Subdomainkopplung."));
    interface->setWhatsThis(interface->toolTip());
    toolbar->addWidget(interface);
    QObject::connect(interface, &QPushButton::clicked, &owner, [this] { addInterfaceCircuit(); });
    auto *remove = new QPushButton(LemEditor::tr("Entfernen"));
    remove->setObjectName(QStringLiteral("lemDelete"));
    remove->setToolTip(
        LemEditor::tr("Entfernt ausgewählte Bauteile oder trennt ausgewählte Leitungen."));
    toolbar->addWidget(remove);
    auto *fit = new QPushButton(LemEditor::tr("Einpassen"));
    fit->setObjectName(QStringLiteral("lemFit"));
    toolbar->addWidget(fit);
    toolbar->addStretch();
    layout->addLayout(toolbar);

    auto *split = new QSplitter;
    auto *palettePanel = new QWidget;
    auto *paletteLayout = new QVBoxLayout(palettePanel);
    paletteLayout->setContentsMargins(0, 0, 0, 0);
    palette = new ComponentPalette;
    palette->setObjectName(QStringLiteral("lemPalette"));
    palette->setToolTip(
        LemEditor::tr("Bauteil auf den Schaltplan ziehen oder auswählen und hinzufügen."));
    palette->setWhatsThis(palette->toolTip());
    for (const auto &entry : paletteEntries) {
        auto *item = new QListWidgetItem(LemEditor::tr(entry.label), palette);
        item->setData(Qt::UserRole, QByteArray(entry.id));
    }
    palette->setCurrentRow(0);
    paletteLayout->addWidget(palette);
    auto *add = new QPushButton(LemEditor::tr("Hinzufügen"));
    add->setObjectName(QStringLiteral("lemAddComponent"));
    paletteLayout->addWidget(add);
    split->addWidget(palettePanel);
    scene = new QGraphicsScene(&owner);
    view = new CircuitView(scene);
    view->setObjectName(QStringLiteral("lemCircuitView"));
    view->setToolTip(LemEditor::tr("Bauteile und Knoten ziehen. Zwei gleichartige Anschlüsse "
                                   "anklicken, um sie zu verbinden. Entf trennt eine ausgewählte "
                                   "Leitung. Strg+Mausrad zoomt von 25 bis 300 Prozent; Einpassen "
                                   "zeigt den gesamten Plan. Anschlüsse sind auch rechts "
                                   "über die Knotenauswahl per Tastatur bearbeitbar."));
    view->setWhatsThis(view->toolTip());
    split->addWidget(view);
    properties = new QScrollArea;
    properties->setObjectName(QStringLiteral("lemProperties"));
    properties->setWidgetResizable(true);
    properties->setMinimumWidth(260);
    split->addWidget(properties);
    split->setStretchFactor(0, 0);
    split->setStretchFactor(1, 1);
    split->setStretchFactor(2, 0);
    split->setSizes({170, 600, 290});
    layout->addWidget(split, 1);
    status = new QLabel;
    status->setObjectName(QStringLiteral("lemStatus"));
    status->setWordWrap(true);
    status->setTextInteractionFlags(Qt::TextSelectableByMouse);
    layout->addWidget(status);
    resultStatus = new QLabel;
    resultStatus->setObjectName(QStringLiteral("lemResultStatus"));
    layout->addWidget(resultStatus);
    results = new QTableWidget(0, 7);
    results->setObjectName(QStringLiteral("lemResults"));
    results->setHorizontalHeaderLabels({LemEditor::tr("Treiber"), "|I| [A RMS]", "|v| [m/s RMS]",
                                        "x [mm Spitze]", "P el. [W]", "P Cu [W]", "Z el. [Ω]"});
    results->setEditTriggers(QAbstractItemView::NoEditTriggers);
    results->setSelectionBehavior(QAbstractItemView::SelectRows);
    results->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    results->horizontalHeader()->setStretchLastSection(true);
    results->setMaximumHeight(150);
    results->setToolTip(LemEditor::tr("Tatsächlich berechnete Treiberdaten der ausgewählten "
                                      "Frequenz. P Cu ist Wärmeleistung, keine Temperatur- "
                                      "oder 10-Sekunden-Belastbarkeitsvorhersage."));
    results->setWhatsThis(results->toolTip());
    layout->addWidget(results);
    portResults = new QTableWidget(0, 4);
    portResults->setObjectName("lemPortResults");
    portResults->setHorizontalHeaderLabels(
        {LemEditor::tr("BEM↔LEM-Anschluss"), "|p| [Pa RMS]", "|Q| [cm³/s RMS]", "P [W]"});
    portResults->setEditTriggers(QAbstractItemView::NoEditTriggers);
    portResults->horizontalHeader()->setSectionResizeMode(QHeaderView::ResizeToContents);
    portResults->horizontalHeader()->setStretchLastSection(true);
    portResults->setMaximumHeight(130);
    portResults->setToolTip(
        LemEditor::tr("Tatsächlich gelöster mittlerer Druck, Volumenstrom und Wirkleistung der "
                      "BEM↔LEM-Anschlüsse. Interface-Q positiv in die Frontluft; keine Kammermoden "
                      "aus dem konzentrierten Volumenmodell."));
    portResults->setWhatsThis(portResults->toolTip());
    layout->addWidget(portResults);

    QObject::connect(enabled, &QCheckBox::toggled, &owner, [this](bool value) {
        if (owner.onEnabledChanged)
            owner.onEnabledChanged(value);
    });
    QObject::connect(preset, &QPushButton::clicked, &owner, [this] { addDriverCircuit(); });
    QObject::connect(add, &QPushButton::clicked, &owner, [this] {
        if (palette->currentItem())
            addComponent(palette->currentItem()->data(Qt::UserRole).toByteArray(),
                         view->mapToScene(view->viewport()->rect().center()));
    });
    QObject::connect(palette, &QListWidget::itemDoubleClicked, &owner,
                     [this](QListWidgetItem *item) {
                         addComponent(item->data(Qt::UserRole).toByteArray(),
                                      view->mapToScene(view->viewport()->rect().center()));
                     });
    QObject::connect(remove, &QPushButton::clicked, &owner, [this] { deleteSelection(); });
    QObject::connect(fit, &QPushButton::clicked, &owner, [this] {
        view->fitInView(scene->itemsBoundingRect().adjusted(-35, -35, 35, 35), Qt::KeepAspectRatio);
    });
    QObject::connect(scene, &QGraphicsScene::selectionChanged, &owner,
                     [this] { selectionChanged(); });
    view->onDrop = [this](auto kind, auto position) { addComponent(kind, position); };
    view->onDelete = [this] { deleteSelection(); };
    view->onCancelWire = [this] {
        pendingTerminal.reset();
        updateStatus();
    };
    rebuildScene();
    updateResults();
}

void LemEditor::Impl::showError(const QString &text) {
    auto colors = status->palette();
    colors.setColor(QPalette::WindowText, QColor(180, 40, 40));
    status->setPalette(colors);
    status->setText(text);
}

void LemEditor::Impl::updateStatus() {
    status->setPalette(owner.palette());
    if (pendingTerminal) {
        status->setText(LemEditor::tr("Zweiten Anschluss oder Knoten wählen · Esc bricht ab."));
        return;
    }
    if (network.elements.empty()) {
        status->setText(LemEditor::tr("Bauteile hinzufügen und Anschlüsse verbinden."));
        return;
    }
    for (const auto &element : network.elements)
        if (const auto *cavity = std::get_if<lem::Cavity>(&element.component);
            cavity && cavity->volume == 0) {
            showError(LemEditor::tr(
                "Luftkammervolumen nicht gesetzt · in cm³ eingeben; noch nicht berechenbar."));
            return;
        }
    try {
        lem::validateNetwork(network);
        for (const auto &part : network.elements) {
            if (const auto *port = std::get_if<lem::BemPort>(&part.component)) {
                for (const auto face : port->faceIds)
                    if (std::none_of(
                            faces.begin(), faces.end(), [face, port](const auto &available) {
                                return available.id == face && portFace(available, port->kind);
                            })) {
                        showError(
                            LemEditor::tr(
                                "BEM-Fläche %1 ist nicht für die gewählte Anschlussart verfügbar "
                                "(Driven oder virtuelles Wall-Interface ohne Apertur).")
                                .arg(face));
                        return;
                    }
            }
        }
        status->setText(LemEditor::tr("Netzwerk gültig · %1 Bauteile · %2 Knoten")
                            .arg(network.elements.size())
                            .arg(network.nodes.size()));
        if (!enabled->isChecked() &&
            std::any_of(network.elements.begin(), network.elements.end(), [](const auto &element) {
                const auto *port = std::get_if<lem::BemPort>(&element.component);
                return port && port->kind == lem::BemPortKind::interfacePort &&
                       !port->faceIds.empty();
            }))
            showError(
                LemEditor::tr("Virtuelles Interface zugeordnet · BEM↔LEM-Kopplung einschalten; die "
                              "Öffnung darf nicht als Wall geschlossen werden."));
    } catch (const std::invalid_argument &error) {
        showError(QString::fromUtf8(error.what()));
    }
}

void LemEditor::Impl::notify(bool rebuild) {
    stale = solution.has_value();
    if (rebuild)
        rebuildScene();
    else
        updateStatus();
    updateResults();
    if (owner.onChanged)
        owner.onChanged();
}

void LemEditor::Impl::rebuildScene() {
    rebuilding = true;
    pendingTerminal.reset();
    scene->clear();
    wires.clear();
    componentItems.clear();
    nodeItems.clear();
    for (std::size_t i = 0; i < network.nodes.size(); ++i) {
        const auto &value = network.nodes[i];
        auto *item = new JunctionItem(i, value.reference, QString::fromStdString(value.name),
                                      domainColor(value.domain));
        item->setToolTip(domainName(value.domain) + QStringLiteral(" · ") +
                         QString::fromStdString(value.name));
        item->setPos(value.position.x, value.position.y);
        scene->addItem(item);
        nodeItems.push_back(item);
        item->onMoving = [this] { updateWires(); };
        item->onMoved = [this](QPointF) { commitPositions(); };
        item->onClicked = [this, id = value.id] {
            if (pendingTerminal) {
                const auto terminal = *pendingTerminal;
                pendingTerminal.reset();
                QTimer::singleShot(0, &owner, [this, terminal, id] {
                    if (const auto *source = terminalNode(terminal))
                        joinNodes(*source, id);
                });
            }
        };
        if (selectedNode && value.id == *selectedNode)
            item->setSelected(true);
    }
    for (std::size_t i = 0; i < network.elements.size(); ++i) {
        auto &part = network.elements[i];
        auto pins = terminals(part);
        const bool driver = std::holds_alternative<lem::Driver>(part.component);
        auto *item = new ComponentItem(i, QString::fromStdString(part.name),
                                       componentName(part.component), driver ? 116 : 86);
        item->setPos(part.position.x, part.position.y);
        scene->addItem(item);
        componentItems.push_back(item);
        item->onMoving = [this] { updateWires(); };
        item->onMoved = [this](QPointF) { commitPositions(); };
        for (std::size_t pin = 0; pin < pins.size(); ++pin) {
            const QPointF point = driver ? QPointF(pin < 2 ? 0 : 200, pin % 2 == 0 ? 70 : 100)
                                         : QPointF(pin == 0 ? 0 : 200, 66);
            const TerminalAddress address{i, pin};
            auto *terminal = new TerminalItem(point, domainColor(terminalDomain(address)),
                                              pins[pin].label, item);
            auto *label = new QGraphicsTextItem(pins[pin].label, item);
            label->setDefaultTextColor(QApplication::palette().color(QPalette::Text));
            label->setPos(point.x() == 0 ? 10 : 105, point.y() - 12);
            terminal->onClicked = [this, address] {
                QTimer::singleShot(0, &owner, [this, address] { terminalClicked(address); });
            };
            auto *wire = new QGraphicsPathItem;
            wire->setFlags(QGraphicsItem::ItemIsSelectable);
            wire->setData(itemKindRole, static_cast<int>(DiagramItemKind::wire));
            wire->setData(itemIndexRole, QVariant::fromValue<qulonglong>(i));
            wire->setData(terminalRole, QVariant::fromValue<qulonglong>(pin));
            wire->setPen(QPen(domainColor(terminalDomain(address)), 2));
            wire->setZValue(-1);
            scene->addItem(wire);
            wires.push_back({address, wire, terminal});
        }
        if (selectedElement && part.id == *selectedElement)
            item->setSelected(true);
    }
    updateWires();
    scene->setSceneRect(scene->itemsBoundingRect().adjusted(-100, -100, 100, 100));
    rebuilding = false;
    updateProperties();
    updateStatus();
}

void LemEditor::Impl::updateWires() {
    for (auto &wire : wires) {
        const auto *id = terminalNode(wire.address);
        if (!id)
            continue;
        const auto it = std::find_if(network.nodes.begin(), network.nodes.end(),
                                     [id](const auto &value) { return value.id == *id; });
        if (it == network.nodes.end()) {
            wire.item->setVisible(false);
            continue;
        }
        const auto index = static_cast<std::size_t>(it - network.nodes.begin());
        const auto start = wire.terminal->scenePos();
        const auto end = nodeItems[index]->scenePos();
        QPainterPath path(start);
        const double middle = (start.x() + end.x()) / 2;
        path.lineTo(middle, start.y());
        path.lineTo(middle, end.y());
        path.lineTo(end);
        wire.item->setPath(path);
        wire.item->setVisible(true);
    }
}

void LemEditor::Impl::commitPositions() {
    // Qt moves every selected item but only the mouse grabber receives release.
    // Commit the complete moved selection as one undoable edit.
    for (std::size_t i = 0; i < network.nodes.size(); ++i)
        network.nodes[i].position = {nodeItems[i]->pos().x(), nodeItems[i]->pos().y()};
    for (std::size_t i = 0; i < network.elements.size(); ++i)
        network.elements[i].position = {componentItems[i]->pos().x(), componentItems[i]->pos().y()};
    notify(false);
}

void LemEditor::Impl::selectionChanged() {
    if (rebuilding)
        return;
    selectedElement.reset();
    selectedNode.reset();
    for (const auto *item : scene->selectedItems()) {
        const auto index = static_cast<std::size_t>(item->data(itemIndexRole).toULongLong());
        switch (static_cast<DiagramItemKind>(item->data(itemKindRole).toInt())) {
        case DiagramItemKind::component:
            if (index < network.elements.size())
                selectedElement = network.elements[index].id;
            break;
        case DiagramItemKind::node:
            if (index < network.nodes.size())
                selectedNode = network.nodes[index].id;
            break;
        default:
            break;
        }
    }
    updateProperties();
    if (selectedElement && owner.onFacesSelected) {
        const auto *selected = element(*selectedElement);
        if (selected)
            if (const auto *port = std::get_if<lem::BemPort>(&selected->component))
                owner.onFacesSelected(port->faceIds);
    }
}

void LemEditor::Impl::terminalClicked(TerminalAddress address) {
    if (!terminalNode(address))
        return;
    if (!pendingTerminal) {
        pendingTerminal = address;
        updateStatus();
        return;
    }
    const auto first = *pendingTerminal;
    pendingTerminal.reset();
    if (first == address) {
        updateStatus();
        return;
    }
    if (terminalDomain(first) != terminalDomain(address)) {
        showError(LemEditor::tr("Anschlüsse verschiedener physikalischer Bereiche lassen sich "
                                "nicht direkt verbinden."));
        return;
    }
    const auto *from = terminalNode(first);
    const auto *to = terminalNode(address);
    if (from && to)
        joinNodes(*from, *to);
}

void LemEditor::Impl::joinNodes(lem::NodeId from, lem::NodeId to) {
    if (from == to)
        return;
    const auto *first = node(from);
    const auto *second = node(to);
    if (!first || !second || first->domain != second->domain) {
        showError(LemEditor::tr("Die Verbindung benötigt zwei Knoten desselben Bereichs."));
        return;
    }
    if (first->reference && !second->reference)
        std::swap(from, to);
    for (auto &part : network.elements)
        for (const auto &pin : terminals(part))
            if (*pin.node == from)
                *pin.node = to;
    std::erase_if(network.nodes, [from](const auto &value) { return value.id == from; });
    if (selectedNode == from)
        selectedNode = to;
    notify();
}

void LemEditor::Impl::attachTerminal(TerminalAddress address, lem::NodeId target) {
    auto *terminal = terminalNode(address);
    const auto *destination = node(target);
    if (!terminal || !destination)
        return;
    if (terminalDomain(address) != destination->domain) {
        showError(LemEditor::tr("Anschluss und Knoten gehören zu verschiedenen Bereichen."));
        return;
    }
    if (*terminal != target) {
        *terminal = target;
        removeUnusedNodes();
        notify();
    }
}

void LemEditor::Impl::removeUnusedNodes() {
    std::set<std::uint32_t> used;
    for (auto &part : network.elements)
        for (const auto &pin : terminals(part))
            used.insert(pin.node->value);
    std::erase_if(network.nodes, [&](const auto &value) {
        return !value.reference && !used.contains(value.id.value);
    });
}

void LemEditor::Impl::addComponent(const QByteArray &kind, QPointF position) {
    const auto domain =
        kind.startsWith("a-") || kind == "cavity" || kind == "bem" || kind == "bem-interface"
            ? lem::Domain::acoustic
        : kind.startsWith("m-") ? lem::Domain::mechanical
                                : lem::Domain::electrical;
    if (kind.endsWith("-ground")) {
        selectedNode = createNode(domain, position, true);
        selectedElement.reset();
        notify();
        return;
    }
    const bool known = std::any_of(paletteEntries.begin(), paletteEntries.end(),
                                   [&kind](const auto &entry) { return kind == entry.id; });
    if (!known)
        return;
    const lem::Terminals pair{createNode(domain, position + QPointF(-65, 66)),
                              createNode(domain, position + QPointF(265, 66))};
    lem::Component component;
    if (kind == "driver") {
        component =
            lem::Driver{pair, createNode(lem::Domain::acoustic, position + QPointF(270, 70)),
                        createNode(lem::Domain::acoustic, position + QPointF(270, 150)),
                        ti100ConceptParameters()};
    } else if (kind == "cavity") {
        component = lem::Cavity{pair, 0};
    } else if (kind == "a-slit") {
        component = lem::AcousticSlit{pair};
    } else if (kind == "a-closed-slit") {
        component = lem::ClosedSlitCavity{pair};
    } else if (kind == "bem" || kind == "bem-interface") {
        component = lem::BemPort{pair,
                                 {},
                                 1,
                                 kind == "bem-interface" ? lem::BemPortKind::interfacePort
                                                         : lem::BemPortKind::drivenSurface};
    } else if (kind.endsWith("-source")) {
        component = lem::EffortSource{pair, {1, 0}};
    } else if (kind.endsWith("-flow")) {
        component = lem::FlowSource{pair, {1, 0}};
    } else if (kind.endsWith("-r")) {
        component = lem::Resistance{pair, domain == lem::Domain::acoustic ? 1000.0 : 1.0};
    } else if (kind.endsWith("-l")) {
        component = lem::Inertance{pair, domain == lem::Domain::acoustic ? 100.0 : .001};
    } else {
        component = lem::Compliance{pair, domain == lem::Domain::electrical ? .000001 : .00001};
    }
    std::uint32_t next = 1;
    for (const auto &part : network.elements)
        next = std::max(next, part.id.value + 1);
    const lem::ElementId id{next};
    const auto name = kind == "driver" ? QStringLiteral("TI 100 · Mmd-Annahme 4 g")
                                       : componentName(component) + QStringLiteral(" %1").arg(next);
    network.elements.push_back(
        {id, name.toStdString(), {position.x(), position.y()}, std::move(component)});
    selectedElement = id;
    selectedNode.reset();
    notify();
}

void LemEditor::Impl::addDriverCircuit() {
    std::set<std::uint32_t> used;
    for (const auto &part : network.elements)
        if (const auto *port = std::get_if<lem::BemPort>(&part.component))
            used.insert(port->faceIds.begin(), port->faceIds.end());
    std::vector<std::uint32_t> assigned;
    for (const auto &face : faces)
        if (face.boundary == BoundaryCondition::driven && !used.contains(face.id)) {
            assigned.push_back(face.id);
            break;
        }
    // The numerical factory requires a face. Keep an explicitly unbound editor port
    // when no unused Driven face is available; validation will prevent simulation.
    auto appended = lem::makeDriverNetwork(ti100ConceptParameters(),
                                           {{assigned.empty() ? 1 : assigned[0]}}, {27, 0}, .00034);
    if (assigned.empty())
        for (auto &part : appended.elements)
            if (auto *port = std::get_if<lem::BemPort>(&part.component))
                port->faceIds.clear();
    const double offsetY = network.elements.empty() ? 0 : scene->itemsBoundingRect().bottom() + 100;
    std::map<std::uint32_t, lem::NodeId> remap;
    for (const auto &value : appended.nodes) {
        if (value.reference) {
            const auto same =
                std::find_if(network.nodes.begin(), network.nodes.end(), [&](const auto &current) {
                    return current.domain == value.domain && current.reference;
                });
            if (same != network.nodes.end()) {
                remap[value.id.value] = same->id;
                continue;
            }
        }
        const auto id = createNode(value.domain, {value.position.x, value.position.y + offsetY},
                                   value.reference);
        remap[value.id.value] = id;
    }
    std::uint32_t next = 1;
    for (const auto &part : network.elements)
        next = std::max(next, part.id.value + 1);
    for (auto &part : appended.elements) {
        part.id = {next++};
        part.name += " " + std::to_string(part.id.value);
        part.position.y += offsetY;
        for (const auto &pin : terminals(part))
            *pin.node = remap.at(pin.node->value);
        if (std::holds_alternative<lem::Driver>(part.component)) {
            part.name = "TI 100 · Mmd-Annahme 4 g";
            selectedElement = part.id;
        }
        network.elements.push_back(std::move(part));
    }
    selectedNode.reset();
    notify();
    view->fitInView(scene->itemsBoundingRect().adjusted(-35, -35, 35, 35), Qt::KeepAspectRatio);
}

void LemEditor::Impl::addInterfaceCircuit() {
    const double y = network.elements.empty() ? 0 : scene->itemsBoundingRect().bottom() + 100;
    const auto pressure = createNode(lem::Domain::acoustic, {260, y + 70});
    const auto reference = createNode(lem::Domain::acoustic, {260, y + 240}, true);
    const lem::Terminals pair{pressure, reference};
    std::set<std::uint32_t> used;
    for (const auto *element : lem::bemPorts(network)) {
        const auto &port = std::get<lem::BemPort>(element->component);
        used.insert(port.faceIds.begin(), port.faceIds.end());
    }
    std::vector<std::uint32_t> assigned;
    for (const auto id : selectedFaces)
        if (!used.contains(id) && std::any_of(faces.begin(), faces.end(), [id](const auto &face) {
                return face.id == id && portFace(face, lem::BemPortKind::interfacePort);
            }))
            assigned.push_back(id);
    std::uint32_t next = 1;
    for (const auto &element : network.elements)
        next = std::max(next, element.id.value + 1);
    network.elements.push_back({{next},
                                "BEM↔LEM-Interface",
                                {0, y},
                                lem::BemPort{pair, assigned, 1, lem::BemPortKind::interfacePort}});
    network.elements.push_back(
        {{next + 1}, "Luftkammer · Volumen offen", {330, y}, lem::Cavity{pair, 0}});
    selectedElement = lem::ElementId{next + 1};
    selectedNode.reset();
    notify();
    if (owner.onFacesSelected)
        owner.onFacesSelected(assigned);
    view->fitInView(scene->itemsBoundingRect().adjusted(-35, -35, 35, 35), Qt::KeepAspectRatio);
}

void LemEditor::Impl::deleteSelection() {
    std::set<std::uint32_t> parts;
    std::set<std::uint32_t> nodes;
    std::vector<TerminalAddress> detached;
    for (const auto *item : scene->selectedItems()) {
        const auto index = static_cast<std::size_t>(item->data(itemIndexRole).toULongLong());
        switch (static_cast<DiagramItemKind>(item->data(itemKindRole).toInt())) {
        case DiagramItemKind::component:
            if (index < network.elements.size())
                parts.insert(network.elements[index].id.value);
            break;
        case DiagramItemKind::node:
            if (index < network.nodes.size())
                nodes.insert(network.nodes[index].id.value);
            break;
        case DiagramItemKind::wire:
            detached.push_back(
                {index, static_cast<std::size_t>(item->data(terminalRole).toULongLong())});
            break;
        }
    }
    if (parts.empty() && nodes.empty() && detached.empty())
        return;
    for (std::size_t i = 0; i < network.elements.size(); ++i) {
        if (parts.contains(network.elements[i].id.value))
            continue;
        auto pins = terminals(network.elements[i]);
        for (std::size_t pin = 0; pin < pins.size(); ++pin) {
            if (nodes.contains(pins[pin].node->value) ||
                std::find(detached.begin(), detached.end(), TerminalAddress{i, pin}) !=
                    detached.end()) {
                const auto domain = terminalDomain({i, pin});
                const auto position = componentItems[i]->pos() + QPointF(-45, 140 + 30 * pin);
                *pins[pin].node = createNode(domain, position);
            }
        }
    }
    std::erase_if(network.elements,
                  [&](const auto &part) { return parts.contains(part.id.value); });
    std::erase_if(network.nodes, [&](const auto &value) { return nodes.contains(value.id.value); });
    removeUnusedNodes();
    selectedElement.reset();
    selectedNode.reset();
    notify();
}

void LemEditor::Impl::updateProperties() {
    if (auto *previous = properties->takeWidget()) {
        for (auto *child : previous->findChildren<QObject *>())
            QObject::disconnect(child, nullptr, &owner, nullptr);
        previous->hide();
        previous->setParent(nullptr);
        previous->deleteLater();
    }
    auto *panel = new QWidget;
    auto *form = new QFormLayout(panel);
    form->setFieldGrowthPolicy(QFormLayout::AllNonFixedFieldsGrow);
    form->setRowWrapPolicy(QFormLayout::WrapLongRows);
    properties->setWidget(panel);
    const auto addNumber = [&](const QString &label, const QString &objectName, double value,
                               double scale, const QString &unit, const QString &description,
                               std::function<void(double)> assign, double minimum = 0.0,
                               double maximum = 1e12, int decimals = 6) {
        auto *input = new DecimalInput;
        input->setObjectName(objectName);
        input->setDecimals(decimals);
        input->setRange(minimum, maximum);
        input->setSingleStep(std::max(std::abs(value * scale) / 20, std::pow(10.0, -decimals)));
        input->setSuffix(unit.isEmpty() ? QString{} : QStringLiteral(" ") + unit);
        input->setValue(value * scale);
        input->setToolTip(description);
        input->setWhatsThis(description);
        form->addRow(label, input);
        QObject::connect(input, &QDoubleSpinBox::valueChanged, &owner,
                         [this, assign = std::move(assign), scale](double displayed) {
                             assign(displayed / scale);
                             notify(false);
                         });
    };
    if (selectedNode && node(*selectedNode)) {
        const auto id = *selectedNode;
        const auto &value = *node(id);
        auto *name = new QLineEdit(QString::fromStdString(value.name));
        name->setObjectName(QStringLiteral("lemNodeName"));
        form->addRow(LemEditor::tr("Knoten"), name);
        form->addRow(LemEditor::tr("Bereich"), new QLabel(domainName(value.domain)));
        auto *reference = new QCheckBox(LemEditor::tr("Referenz (Potential = 0)"));
        reference->setObjectName(QStringLiteral("lemNodeReference"));
        reference->setChecked(value.reference);
        reference->setToolTip(
            LemEditor::tr("Expliziter Nullpunkt dieses physikalischen Bereichs. "
                          "Jede zusammenhängende Gruppe benötigt eine Referenz."));
        reference->setWhatsThis(reference->toolTip());
        form->addRow(reference);
        QObject::connect(name, &QLineEdit::editingFinished, &owner, [this, id, name] {
            if (auto *current = node(id); current && current->name != name->text().toStdString()) {
                current->name = name->text().toStdString();
                notify();
            }
        });
        QObject::connect(reference, &QCheckBox::toggled, &owner, [this, id](bool checked) {
            if (auto *current = node(id)) {
                current->reference = checked;
                notify();
            }
        });
        return;
    }
    if (!selectedElement || !element(*selectedElement)) {
        auto *hint = new QLabel(LemEditor::tr("Bauteil oder Knoten auswählen."));
        hint->setWordWrap(true);
        form->addRow(hint);
        return;
    }
    const auto id = *selectedElement;
    auto &part = *element(id);
    const auto elementIndex = static_cast<std::size_t>(&part - network.elements.data());
    auto *name = new QLineEdit(QString::fromStdString(part.name));
    name->setObjectName(QStringLiteral("lemElementName"));
    form->addRow(LemEditor::tr("Name"), name);
    form->addRow(LemEditor::tr("Bauteil"), new QLabel(componentName(part.component)));
    QObject::connect(name, &QLineEdit::editingFinished, &owner, [this, id, name] {
        if (auto *current = element(id); current && current->name != name->text().toStdString()) {
            current->name = name->text().toStdString();
            notify();
        }
    });
    auto pins = terminals(part);
    for (std::size_t pin = 0; pin < pins.size(); ++pin) {
        const auto domain = terminalDomain({elementIndex, pin});
        auto *connection = new QComboBox;
        connection->setObjectName(QStringLiteral("lemTerminal%1").arg(pin));
        connection->addItem(LemEditor::tr("Neuer Knoten"), 0U);
        for (const auto &value : network.nodes)
            if (value.domain == domain)
                connection->addItem(QString::fromStdString(value.name) +
                                        (value.reference ? QStringLiteral(" [0]") : QString{}),
                                    value.id.value);
        connection->setCurrentIndex(connection->findData(pins[pin].node->value));
        connection->setToolTip(LemEditor::tr("%1: Anschluss an einen Knoten desselben Bereichs. "
                                             "Ein neuer Knoten trennt diesen Anschluss ab.")
                                   .arg(domainName(domain)));
        connection->setWhatsThis(connection->toolTip());
        form->addRow(pins[pin].label, connection);
        QObject::connect(
            connection, &QComboBox::activated, &owner,
            [this, id, pin, connection, domain](int index) {
                auto *current = element(id);
                if (!current)
                    return;
                auto target = lem::NodeId{connection->itemData(index).toUInt()};
                if (!target.value)
                    target = createNode(
                        domain, QPointF(current->position.x - 60,
                                        current->position.y + 140 + 30 * static_cast<double>(pin)));
                attachTerminal({static_cast<std::size_t>(current - network.elements.data()), pin},
                               target);
            });
    }
    const auto domain = terminalDomain({elementIndex, 0});
    std::visit(
        [&](const auto &component) {
            using Part = std::decay_t<decltype(component)>;
            if constexpr (std::is_same_v<Part, lem::Driver>) {
                const auto parameter = [&](const QString &label, const char *objectName,
                                           double lem::DriverParameters::*member, double scale,
                                           const QString &unit, const QString &description) {
                    addNumber(label, QString::fromLatin1(objectName), component.parameters.*member,
                              scale, unit, description, [this, id, member](double value) {
                                  if (auto *current = element(id))
                                      std::get<lem::Driver>(current->component).parameters.*member =
                                          value;
                              });
                };
                parameter("Re", "lemDriverResistance", &lem::DriverParameters::resistance, 1, "Ω",
                          LemEditor::tr("Gleichstromwiderstand der Schwingspule [Ω]."));
                parameter("Le", "lemDriverInductance", &lem::DriverParameters::inductance, 1000,
                          "mH",
                          LemEditor::tr("Konstante Schwingspuleninduktivität [mH], intern H."));
                parameter("Bl", "lemDriverForceFactor", &lem::DriverParameters::forceFactor, 1,
                          "N/A",
                          LemEditor::tr("Kraftfaktor des elektrodynamischen Wandlers [N/A]."));
                parameter(
                    "Mmd", "lemDriverMovingMass", &lem::DriverParameters::movingMass, 1000, "g",
                    LemEditor::tr("Trockene bewegte Masse ohne externe Luftlast [g], intern kg. "
                                  "Die volle BEM-Last liefert die äußere Luftlast. "
                                  "Im TI-100-Startkreis sind 4 g aus Mms nur eine Annahme; "
                                  "Mmd muss verifiziert oder ersetzt werden."));
                parameter("Cms", "lemDriverCompliance", &lem::DriverParameters::compliance, 1000,
                          "mm/N", LemEditor::tr("Mechanische Nachgiebigkeit [mm/N], intern m/N."));
                parameter("Rms", "lemDriverMechanicalResistance",
                          &lem::DriverParameters::mechanicalResistance, 1, "N·s/m",
                          LemEditor::tr("Mechanischer Verlustwiderstand [N·s/m]."));
                parameter("Sd", "lemDriverPistonArea", &lem::DriverParameters::pistonArea, 10000,
                          "cm²",
                          LemEditor::tr("Effektive Membranfläche [cm²], intern m². "
                                        "Volumenstrom = Sd × Membrangeschwindigkeit."));
            } else if constexpr (std::is_same_v<Part, lem::AcousticSlit>) {
                const auto parameter = [&](const QString &label, const QString &name,
                                           double lem::AcousticSlit::*member, double scale,
                                           const QString &unit, const QString &description) {
                    addNumber(label, name, component.*member, scale, unit, description,
                              [this, id, member](double value) {
                                  if (auto *current = element(id))
                                      std::get<lem::AcousticSlit>(current->component).*member =
                                          value;
                              });
                };
                parameter(LemEditor::tr("Länge"), "lemSlitLength", &lem::AcousticSlit::length, 1000,
                          "mm",
                          LemEditor::tr("Reale Kanallänge zwischen den zwei akustischen Knoten "
                                        "[mm], intern m. Der Spaltkanal bildet den LEM-Anteil "
                                        "hinter dem mathematischen BEM-Interface ab."));
                parameter(
                    LemEditor::tr("Querschnitt"), "lemSlitArea", &lem::AcousticSlit::area, 1e6,
                    "mm²",
                    LemEditor::tr("Freier Kanalquerschnitt senkrecht zur Strömung [mm²], intern "
                                  "m². Das BEM-Interface ist die Kopplungsgrenze; bereits in BEM "
                                  "aufgelöste Kanalabschnitte hier nicht erneut zählen."));
                parameter(
                    LemEditor::tr("Spaltweite"), "lemSlitGap", &lem::AcousticSlit::gap, 1000, "mm",
                    LemEditor::tr(
                        "Volle lichte Weite des ebenen Spalts [mm], intern m; keine Halbweite. Sie "
                        "bestimmt das dynamische viskose Geschwindigkeitsprofil."));
                parameter(
                    LemEditor::tr("Viskosität μ"), "lemSlitViscosity",
                    &lem::AcousticSlit::dynamicViscosity, 1e6, "µPa·s",
                    LemEditor::tr("Dynamische Viskosität des Mediums [µPa·s], intern Pa·s. Der "
                                  "Zweipol berechnet die frequenzabhängige viskose Spaltimpedanz. "
                                  "Alle vier Parameter müssen für die Berechnung positiv sein."));
                for (const auto *name :
                     {"lemSlitLength", "lemSlitArea", "lemSlitGap", "lemSlitViscosity"})
                    panel->findChild<QDoubleSpinBox *>(name)->setSpecialValueText(
                        LemEditor::tr("nicht gesetzt"));
            } else if constexpr (std::is_same_v<Part, lem::ClosedSlitCavity>) {
                auto *volume = new QLabel;
                volume->setObjectName("lemClosedSlitVolume");
                volume->setToolTip(
                    LemEditor::tr("Geometrisches Volumen A × L der geschlossenen 1D-Spaltleitung "
                                  "[cm³ = mL]. Aus Länge und Querschnitt berechnet; kein "
                                  "zusätzliches konzentriertes Kammervolumen."));
                volume->setWhatsThis(volume->toolTip());
                const auto updateVolume = [this, id, volume] {
                    if (const auto *current = element(id)) {
                        const auto &value = std::get<lem::ClosedSlitCavity>(current->component);
                        const auto cubicCentimetres = value.length * value.area * 1e6;
                        volume->setText(!std::isfinite(cubicCentimetres) ? LemEditor::tr("ungültig")
                                        : cubicCentimetres > 0
                                            ? QString::number(cubicCentimetres, 'g', 8) + " cm³"
                                            : LemEditor::tr("nicht gesetzt"));
                    }
                };
                const auto parameter = [&](const QString &label, const QString &name,
                                           double lem::ClosedSlitCavity::*member, double scale,
                                           const QString &unit, const QString &description,
                                           double minimum = 0.) {
                    addNumber(
                        label, name, component.*member, scale, unit, description,
                        [this, id, member, updateVolume](double value) {
                            if (auto *current = element(id))
                                std::get<lem::ClosedSlitCavity>(current->component).*member = value;
                            updateVolume();
                        },
                        minimum);
                };
                parameter(LemEditor::tr("Länge"), "lemClosedSlitLength",
                          &lem::ClosedSlitCavity::length, 1000, "mm",
                          LemEditor::tr("Axiale Tiefe der am Ende geschlossenen, gleichförmigen "
                                        "Spaltleitung [mm], intern m. Das 1D-Modell berücksichtigt "
                                        "viskose und thermische Wellen entlang dieser Länge."));
                parameter(
                    LemEditor::tr("Querschnitt"), "lemClosedSlitArea", &lem::ClosedSlitCavity::area,
                    1e6, "mm²",
                    LemEditor::tr("Freier Querschnitt der geschlossenen Spaltleitung senkrecht zur "
                                  "Strömung [mm²], intern m². Ihr Volumen ist A × L. Bereits in "
                                  "BEM aufgelöste Abschnitte hier nicht erneut zählen."));
                parameter(LemEditor::tr("Volle Breite"), "lemClosedSlitGap",
                          &lem::ClosedSlitCavity::gap, 1000, "mm",
                          LemEditor::tr("Volle lichte Weite zwischen den beiden Spaltwänden [mm], "
                                        "intern m; keine Halbweite. Sie bestimmt die viskose und "
                                        "thermische Grenzschicht der 1D-Welle."));
                parameter(LemEditor::tr("Viskosität μ"), "lemClosedSlitViscosity",
                          &lem::ClosedSlitCavity::dynamicViscosity, 1e6, "µPa·s",
                          LemEditor::tr(
                              "Positive dynamische Viskosität des Mediums [µPa·s], intern Pa·s."));
                parameter(
                    LemEditor::tr("Wärmekapazität γ"), "lemClosedSlitGamma",
                    &lem::ClosedSlitCavity::heatCapacityRatio, 1, "[-]",
                    LemEditor::tr("Dimensionsloses Wärmekapazitätsverhältnis cp/cv des Mediums, "
                                  "mindestens 1. Es bestimmt mit der Prandtl-Zahl die thermische "
                                  "Kompressibilität der geschlossenen 1D-Spaltleitung."),
                    1);
                parameter(LemEditor::tr("Prandtl-Zahl"), "lemClosedSlitPrandtl",
                          &lem::ClosedSlitCavity::prandtlNumber, 1, "[-]",
                          LemEditor::tr(
                              "Dimensionslose positive Prandtl-Zahl des Mediums. Zusammen mit γ "
                              "und μ bestimmt sie die thermische Grenzschicht im 1D-Modell."));
                updateVolume();
                form->addRow(LemEditor::tr("Volumen A·L"), volume);
                for (const auto *name :
                     {"lemClosedSlitLength", "lemClosedSlitArea", "lemClosedSlitGap",
                      "lemClosedSlitViscosity", "lemClosedSlitPrandtl"})
                    panel->findChild<QDoubleSpinBox *>(name)->setSpecialValueText(
                        LemEditor::tr("nicht gesetzt"));
            } else if constexpr (std::is_same_v<Part, lem::Cavity>) {
                addNumber(LemEditor::tr("Volumen"), "lemVolume", component.volume, 1e6, "cm³",
                          LemEditor::tr(
                              "Geschlossenes, konzentriertes Luftvolumen [cm³ = mL]. Intern m³. "
                              "Nachgiebigkeit V/(ρc²); keine verteilten Kammermoden. 0 bedeutet "
                              "ungesetzt und sperrt die Berechnung. Ein als Annahme benannter Wert "
                              "ist zu verifizieren."),
                          [this, id](double value) {
                              if (auto *current = element(id))
                                  std::get<lem::Cavity>(current->component).volume = value;
                          });
                panel->findChild<QDoubleSpinBox *>("lemVolume")
                    ->setSpecialValueText(LemEditor::tr("nicht gesetzt"));
            } else if constexpr (std::is_same_v<Part, lem::BemPort>) {
                if (component.needsAssignmentReview) {
                    auto *review = new QLabel(
                        LemEditor::tr("Zuordnung nach Geometrieänderung prüfen · passende Flächen "
                                      "neu übernehmen oder Anschluss entfernen."));
                    review->setObjectName("lemPortAssignmentReview");
                    review->setWordWrap(true);
                    review->setStyleSheet("color: #a72e2b");
                    form->addRow(review);
                }
                auto *kind = new QComboBox;
                kind->setObjectName("lemPortKind");
                kind->addItem(LemEditor::tr("Antriebsfläche (Driven)"));
                kind->addItem(LemEditor::tr("Virtuelles BEM↔LEM-Interface"));
                kind->setCurrentIndex(component.kind == lem::BemPortKind::interfacePort ? 1 : 0);
                kind->setToolTip(LemEditor::tr(
                    "Interface: mathematischer Übergang von der BEM-Frontluft zum "
                    "akustischen LEM-Netz. Geometrisch Wall, physikalisch offener Anschluss; keine "
                    "Materialwand und keine FEM-/BEM↔BEM-Kopplung. Wechsel der Art löst die "
                    "bisherige Flächenzuordnung."));
                kind->setWhatsThis(kind->toolTip());
                form->addRow(LemEditor::tr("Anschlussart"), kind);
                QObject::connect(kind, &QComboBox::currentIndexChanged, &owner,
                                 [this, id](int index) {
                                     if (auto *current = element(id)) {
                                         auto &port = std::get<lem::BemPort>(current->component);
                                         port.kind = index == 1 ? lem::BemPortKind::interfacePort
                                                                : lem::BemPortKind::drivenSurface;
                                         port.needsAssignmentReview |= !port.faceIds.empty();
                                         port.faceIds.clear();
                                         port.velocityWeights.clear();
                                         notify();
                                     }
                                 });
                auto *orientation = new QComboBox;
                orientation->setObjectName(QStringLiteral("lemPortPolarity"));
                orientation->addItem(component.kind == lem::BemPortKind::interfacePort
                                         ? LemEditor::tr("+ in BEM-Frontluft")
                                         : LemEditor::tr("+ aus Material"),
                                     1);
                orientation->addItem(component.kind == lem::BemPortKind::interfacePort
                                         ? LemEditor::tr("− aus BEM-Frontluft")
                                         : LemEditor::tr("− in Material"),
                                     -1);
                orientation->setCurrentIndex(component.polarity == 1 ? 0 : 1);
                orientation->setToolTip(LemEditor::tr(
                    "Q BEM = Vorzeichen × Q Netzwerk. "
                    "Interface: positives Q injiziert in die BEM-Frontluft, entgegen der aus AIR "
                    "gerichteten Netznormale. Antriebsport: positiv aus dem Material in die Luft. "
                    "Der negative Anschluss muss Druckreferenz sein."));
                orientation->setWhatsThis(orientation->toolTip());
                form->addRow(LemEditor::tr("Orientierung"), orientation);
                QObject::connect(orientation, &QComboBox::currentIndexChanged, &owner,
                                 [this, id, orientation](int) {
                                     if (auto *current = element(id)) {
                                         std::get<lem::BemPort>(current->component).polarity =
                                             orientation->currentData().toInt();
                                         notify(false);
                                     }
                                 });
                auto *selection = new QListWidget;
                selection->setObjectName(QStringLiteral("lemPortFaces"));
                selection->setMinimumHeight(110);
                selection->setMaximumHeight(220);
                selection->setToolTip(LemEditor::tr(
                    "Eine oder mehrere zur Anschlussart passende Flächen bilden einen "
                    "gemeinsamen BEM-Port. Interfaceflächen bewegen Luft gleichmäßig; "
                    "Driven-Flächen verwenden ihre relativen Geschwindigkeitsgewichte. Jede Fläche "
                    "darf genau "
                    "einem Port zugeordnet werden."));
                selection->setWhatsThis(selection->toolTip());
                const auto addFace = [&](std::uint32_t face, const QString &label, bool checked) {
                    auto *item = new QListWidgetItem(label, selection);
                    item->setData(Qt::UserRole, face);
                    item->setFlags(item->flags() | Qt::ItemIsUserCheckable);
                    item->setCheckState(checked ? Qt::Checked : Qt::Unchecked);
                };
                for (const auto &face : faces)
                    if (portFace(face, component.kind))
                        addFace(face.id, face.label,
                                std::find(component.faceIds.begin(), component.faceIds.end(),
                                          face.id) != component.faceIds.end());
                for (const auto face : component.faceIds)
                    if (std::none_of(faces.begin(), faces.end(), [&](const auto &item) {
                            return item.id == face && portFace(item, component.kind);
                        }))
                        addFace(face, LemEditor::tr("Fläche %1 (nicht verfügbar)").arg(face), true);
                form->addRow(LemEditor::tr("BEM-Flächen"), selection);
                if (component.kind == lem::BemPortKind::drivenSurface) {
                    auto *weights = new QLineEdit(portWeightsText(component.velocityWeights));
                    weights->setObjectName("lemPortVelocityWeights");
                    weights->setPlaceholderText(LemEditor::tr("ID:Faktor; … · leer = 1"));
                    weights->setToolTip(LemEditor::tr(
                        "Relative Normalgeschwindigkeit je zugeordneter Driven-Fläche, "
                        "dimensionslos und positiv. Format ID:Faktor; ID:Faktor mit Dezimalpunkt. "
                        "Fehlende IDs haben Faktor 1. v_f = w_f Q / Σ(w A); der Portdruck wird mit "
                        "denselben Gewichten gemittelt. Beispiel 12:0.5 = halbe Geschwindigkeit "
                        "dieser Fläche bei gleicher Phase. Interfaceports bleiben gleichförmig."));
                    weights->setWhatsThis(weights->toolTip());
                    form->addRow(LemEditor::tr("vₙ-Gewichte"), weights);
                    auto *apply = new QPushButton(LemEditor::tr("Gewichte übernehmen"));
                    apply->setObjectName("lemApplyVelocityWeights");
                    form->addRow(apply);
                    auto *error = new QLabel;
                    error->setObjectName("lemPortWeightsError");
                    error->setWordWrap(true);
                    error->setStyleSheet("color: #a72e2b");
                    form->addRow(error);
                    const auto commit = [this, id, weights, error] {
                        auto *current = element(id);
                        if (!current)
                            return;
                        auto &port = std::get<lem::BemPort>(current->component);
                        try {
                            auto parsed = parsedPortWeights(weights->text(), port);
                            error->clear();
                            weights->setText(portWeightsText(parsed));
                            if (port.velocityWeights == parsed)
                                return;
                            port.velocityWeights = std::move(parsed);
                            notify(false);
                        } catch (const std::invalid_argument &failure) {
                            error->setText(QString::fromUtf8(failure.what()));
                        }
                    };
                    QObject::connect(apply, &QPushButton::clicked, &owner, commit);
                    QObject::connect(weights, &QLineEdit::returnPressed, &owner, commit);
                }
                auto *takeSelection = new QPushButton(LemEditor::tr("3D-Auswahl übernehmen"));
                takeSelection->setObjectName("lemTakeSelectedFaces");
                takeSelection->setToolTip(
                    LemEditor::tr("Übernimmt passende Flächen aus der aktuellen 3D-/CAD-Auswahl "
                                  "für diesen Port. Apertur und bereits anderweitig zugeordnete "
                                  "Flächen werden nicht übernommen."));
                takeSelection->setWhatsThis(takeSelection->toolTip());
                takeSelection->setEnabled(!selectedFaces.empty());
                form->addRow(takeSelection);
                QObject::connect(takeSelection, &QPushButton::clicked, &owner, [this, id] {
                    auto *current = element(id);
                    if (!current)
                        return;
                    auto &port = std::get<lem::BemPort>(current->component);
                    const bool wasBound = !port.faceIds.empty();
                    std::set<std::uint32_t> used;
                    for (const auto *other : lem::bemPorts(network))
                        if (other->id != id) {
                            const auto &mapped = std::get<lem::BemPort>(other->component).faceIds;
                            used.insert(mapped.begin(), mapped.end());
                        }
                    port.faceIds.clear();
                    for (const auto faceId : selectedFaces)
                        if (!used.contains(faceId) &&
                            std::any_of(faces.begin(), faces.end(), [&](const auto &face) {
                                return face.id == faceId && portFace(face, port.kind);
                            }))
                            port.faceIds.push_back(faceId);
                    retainPortWeights(port);
                    if (!port.faceIds.empty())
                        port.needsAssignmentReview = false;
                    else
                        port.needsAssignmentReview |=
                            wasBound && port.kind == lem::BemPortKind::interfacePort;
                    if (owner.onFacesSelected)
                        owner.onFacesSelected(port.faceIds);
                    notify();
                });
                QObject::connect(
                    selection, &QListWidget::itemChanged, &owner,
                    [this, id, selection](QListWidgetItem *) {
                        if (auto *current = element(id)) {
                            auto &port = std::get<lem::BemPort>(current->component);
                            const bool wasBound = !port.faceIds.empty();
                            auto &assigned = port.faceIds;
                            assigned.clear();
                            for (int i = 0; i < selection->count(); ++i)
                                if (selection->item(i)->checkState() == Qt::Checked)
                                    assigned.push_back(
                                        selection->item(i)->data(Qt::UserRole).toUInt());
                            retainPortWeights(port);
                            const bool validAssignment =
                                !assigned.empty() &&
                                std::all_of(
                                    assigned.begin(), assigned.end(), [&](const auto faceId) {
                                        const bool available = std::any_of(
                                            faces.begin(), faces.end(), [&](const auto &face) {
                                                return face.id == faceId &&
                                                       portFace(face, port.kind);
                                            });
                                        const bool used = std::any_of(
                                            network.elements.begin(), network.elements.end(),
                                            [&](const auto &other) {
                                                const auto *otherPort =
                                                    std::get_if<lem::BemPort>(&other.component);
                                                return other.id != id && otherPort &&
                                                       std::find(otherPort->faceIds.begin(),
                                                                 otherPort->faceIds.end(),
                                                                 faceId) !=
                                                           otherPort->faceIds.end();
                                            });
                                        return available && !used;
                                    });
                            if (validAssignment)
                                port.needsAssignmentReview = false;
                            else
                                port.needsAssignmentReview |=
                                    wasBound && port.kind == lem::BemPortKind::interfacePort;
                            if (owner.onFacesSelected)
                                owner.onFacesSelected(assigned);
                            notify(false);
                        }
                    });
            } else if constexpr (std::is_same_v<Part, lem::EffortSource> ||
                                 std::is_same_v<Part, lem::FlowSource>) {
                constexpr bool effort = std::is_same_v<Part, lem::EffortSource>;
                const QString unit =
                    domain == lem::Domain::electrical   ? (effort ? "V RMS" : "A RMS")
                    : domain == lem::Domain::mechanical ? (effort ? "N RMS" : "m/s RMS")
                                                        : (effort ? "Pa RMS" : "m³/s RMS");
                addNumber(LemEditor::tr("Amplitude"), "lemAmplitude", std::abs(component.amplitude),
                          1, unit,
                          LemEditor::tr("Effektivwert der harmonischen Quelle. "
                                        "Zeitkonvention exp(+jωt), positive Richtung + → −."),
                          [this, id](double value) {
                              if (auto *current = element(id)) {
                                  auto &amplitude = std::get<Part>(current->component).amplitude;
                                  amplitude = std::polar(value, std::arg(amplitude));
                              }
                          });
                addNumber(
                    LemEditor::tr("Phase"), "lemPhase", std::arg(component.amplitude),
                    180 / std::numbers::pi, "°",
                    LemEditor::tr("Phase der RMS-Quelle in Grad für exp(+jωt)."),
                    [this, id](double value) {
                        if (auto *current = element(id)) {
                            auto &amplitude = std::get<Part>(current->component).amplitude;
                            amplitude = std::polar(std::abs(amplitude), value);
                        }
                    },
                    -180, 180, 3);
            } else {
                const bool electric = domain == lem::Domain::electrical;
                const bool mechanical = domain == lem::Domain::mechanical;
                if constexpr (std::is_same_v<Part, lem::Resistance>) {
                    addNumber("R", "lemResistance", component.resistance, 1,
                              electric     ? "Ω"
                              : mechanical ? "N·s/m"
                                           : "Pa·s/m³",
                              LemEditor::tr(
                                  "Positiver Verlustwiderstand; Potentialdifferenz = R × Fluss."),
                              [this, id](double value) {
                                  if (auto *current = element(id))
                                      std::get<Part>(current->component).resistance = value;
                              });
                } else if constexpr (std::is_same_v<Part, lem::Inertance>) {
                    addNumber(electric ? "L" : "M", "lemInertance", component.inertance,
                              electric || mechanical ? 1000 : 1,
                              electric     ? "mH"
                              : mechanical ? "g"
                                           : "Pa·s²/m³",
                              LemEditor::tr("Induktivität bzw. träge Masse: Z = jωM. "
                                            "Elektrisch intern H, mechanisch kg."),
                              [this, id](double value) {
                                  if (auto *current = element(id))
                                      std::get<Part>(current->component).inertance = value;
                              });
                } else if constexpr (std::is_same_v<Part, lem::Compliance>) {
                    addNumber("C", "lemCompliance", component.compliance,
                              electric     ? 1e6
                              : mechanical ? 1000
                                           : 1e9,
                              electric     ? "µF"
                              : mechanical ? "mm/N"
                                           : "mm³/Pa",
                              LemEditor::tr("Kapazität bzw. Nachgiebigkeit: Z = 1/(jωC). "
                                            "Intern F, m/N oder m³/Pa."),
                              [this, id](double value) {
                                  if (auto *current = element(id))
                                      std::get<Part>(current->component).compliance = value;
                              });
                }
            }
        },
        part.component);
}

void LemEditor::Impl::updateResults() {
    portResults->setRowCount(solution ? static_cast<int>(solution->ports.size()) : 0);
    portResults->setVisible(solution && !solution->ports.empty());
    if (solution)
        for (std::size_t row = 0; row < solution->ports.size(); ++row) {
            const auto &port = solution->ports[row];
            const auto *part = element(port.id);
            const std::array<QString, 4> cells{
                part ? QString::fromStdString(part->name) : QString::number(port.id.value),
                QString::number(std::abs(port.pressure), 'g', 6),
                QString::number(std::abs(port.volumeVelocity) * 1e6, 'g', 6),
                QString::number(port.activePower, 'g', 6)};
            for (std::size_t column = 0; column < cells.size(); ++column) {
                auto *item = new QTableWidgetItem(cells[column]);
                if (stale)
                    item->setForeground(owner.palette().color(QPalette::Disabled, QPalette::Text));
                portResults->setItem(static_cast<int>(row), static_cast<int>(column), item);
            }
        }
    results->setRowCount(solution ? static_cast<int>(solution->drivers.size()) : 0);
    if (!solution) {
        resultStatus->setText(LemEditor::tr("Keine berechneten Treiberdaten."));
        return;
    }
    resultStatus->setText((stale ? LemEditor::tr("Ergebnisse veraltet · %1 Hz · Residuum %2")
                                 : LemEditor::tr("BEM↔LEM-Ergebnisse · %1 Hz · Residuum %2"))
                              .arg(solution->frequency, 0, 'g', 7)
                              .arg(solution->relativeResidual, 0, 'g', 3));
    for (std::size_t row = 0; row < solution->drivers.size(); ++row) {
        const auto &driver = solution->drivers[row];
        const auto *part = element(driver.id);
        const auto electrical =
            std::find_if(solution->elements.begin(), solution->elements.end(),
                         [&](const auto &value) { return value.id == driver.id; });
        const auto impedance =
            electrical != solution->elements.end() && std::abs(driver.current) > 0
                ? complexText(electrical->effort / driver.current)
                : QStringLiteral("—");
        const std::array<QString, 7> cells{
            part ? QString::fromStdString(part->name)
                 : LemEditor::tr("Treiber %1").arg(driver.id.value),
            QString::number(std::abs(driver.current), 'g', 5),
            QString::number(std::abs(driver.velocity), 'g', 5),
            QString::number(std::sqrt(2.0) * std::abs(driver.displacement) * 1000, 'g', 5),
            QString::number(driver.electricalPower, 'g', 5),
            QString::number(driver.copperPower, 'g', 5),
            impedance};
        for (std::size_t column = 0; column < cells.size(); ++column) {
            auto *item = new QTableWidgetItem(cells[column]);
            if (stale)
                item->setForeground(owner.palette().color(QPalette::Disabled, QPalette::Text));
            results->setItem(static_cast<int>(row), static_cast<int>(column), item);
        }
    }
}

LemEditor::LemEditor(QWidget *parent) : QWidget(parent), impl_(std::make_unique<Impl>(*this)) {
    setObjectName(QStringLiteral("lemEditor"));
    impl_->createUi();
}

LemEditor::~LemEditor() {
    for (auto *child : findChildren<QObject *>())
        QObject::disconnect(child, nullptr, this, nullptr);
}

const lem::Network &LemEditor::network() const {
    return impl_->network;
}

void LemEditor::setNetwork(lem::Network network) {
    impl_->network = std::move(network);
    impl_->selectedElement.reset();
    impl_->selectedNode.reset();
    impl_->rebuildScene();
}

void LemEditor::setFaces(std::vector<LemFaceItem> faces) {
    if (impl_->faces == faces)
        return;
    impl_->faces = std::move(faces);
    impl_->updateProperties();
    impl_->updateStatus();
}

void LemEditor::setSelectedFaces(std::vector<std::uint32_t> faces) {
    if (impl_->selectedFaces == faces)
        return;
    impl_->selectedFaces = std::move(faces);
    impl_->updateProperties();
}

bool LemEditor::couplingEnabled() const {
    return impl_->enabled->isChecked();
}

void LemEditor::setCouplingEnabled(bool enabled) {
    const QSignalBlocker blocker(impl_->enabled);
    impl_->enabled->setChecked(enabled);
}

void LemEditor::setSolution(std::optional<lem::Solution> solution, bool stale) {
    impl_->solution = std::move(solution);
    impl_->stale = stale;
    impl_->updateResults();
}

} // namespace hornsculpt
