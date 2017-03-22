#include "godel_simple_gui/states/simulating_state.h"
// prev
#include "godel_simple_gui/states/wait_to_simulate_state.h"
// next
#include "godel_simple_gui/states/wait_to_execute_state.h"

#include <ros/console.h>
#include "godel_simple_gui/blending_widget.h"

#include "godel_msgs/SelectMotionPlan.h"

#include <QtConcurrent/QtConcurrentRun>

const static std::string SELECT_MOTION_PLAN_ACTION_SERVER_NAME = "select_motion_plan_as";

godel_simple_gui::SimulatingState::SimulatingState(const std::vector<std::string>& plans)
    : plan_names_(plans), select_motion_plan_action_client_(SELECT_MOTION_PLAN_ACTION_SERVER_NAME, true)
{
}

void godel_simple_gui::SimulatingState::onStart(BlendingWidget& gui)
{
  gui.setText("Simulating...");
  gui.setButtonsEnabled(false);

  QtConcurrent::run(this, &SimulatingState::simulateAll);
}

void godel_simple_gui::SimulatingState::onExit(BlendingWidget& gui) { gui.setButtonsEnabled(true); }

// Handlers for the fixed buttons
void godel_simple_gui::SimulatingState::onNext(BlendingWidget& gui) {}

void godel_simple_gui::SimulatingState::onBack(BlendingWidget& gui) {}

void godel_simple_gui::SimulatingState::onReset(BlendingWidget& gui) {}

void godel_simple_gui::SimulatingState::simulateAll()
{
  for (std::size_t i = 0; i < plan_names_.size(); ++i)
  {
    simulate_next_ = false;
    simulateOne(plan_names_[i]);
    while(!simulate_next_);  // wait for the last simulation to finish to avoid pre-empting the plan
  }

  Q_EMIT newStateAvailable(new WaitToExecuteState(plan_names_));
}

void godel_simple_gui::SimulatingState::simulateOne(const std::string& plan)
{
  godel_msgs::SelectMotionPlanActionGoal goal;
  goal.goal.name = plan;
  goal.goal.simulate = true;
  goal.goal.wait_for_execution = true;
  select_motion_plan_action_client_.sendGoal(goal.goal,
    boost::bind(&godel_simple_gui::SimulatingState::selectMotionPlanDoneCallback, this, _1, _2));
}

void godel_simple_gui::SimulatingState::selectMotionPlanDoneCallback(const actionlib::SimpleClientGoalState& state,
            const godel_msgs::SelectMotionPlanResultConstPtr& result)
{
  simulate_next_ = true;
}
