#include "godel_simple_gui/states/executing_state.h"
// prev
#include "godel_simple_gui/states/wait_to_execute_state.h"
// There is no next!
// error
#include "godel_simple_gui/states/error_state.h"

#include <ros/console.h>
#include "godel_simple_gui/blending_widget.h"

#include "godel_msgs/SelectMotionPlan.h"

#include <QtConcurrent/QtConcurrentRun>

const static std::string SELECT_MOTION_PLAN_ACTION_SERVER_NAME = "select_motion_plan_as";

struct BadExecutionError
{
  BadExecutionError(const std::string& what) : what(what) {}
  std::string what;
};

godel_simple_gui::ExecutingState::ExecutingState(const std::vector<std::string>& plans)
    : plan_names_(plans), select_motion_plan_action_client_(SELECT_MOTION_PLAN_ACTION_SERVER_NAME, true)
{
}

void godel_simple_gui::ExecutingState::onStart(BlendingWidget& gui)
{
  gui.setText("Executing...");
  gui.setButtonsEnabled(false);

  QtConcurrent::run(this, &ExecutingState::executeAll);
}

void godel_simple_gui::ExecutingState::onExit(BlendingWidget& gui) { gui.setButtonsEnabled(true); }

// Handlers for the fixed buttons
void godel_simple_gui::ExecutingState::onNext(BlendingWidget& gui) {}

void godel_simple_gui::ExecutingState::onBack(BlendingWidget& gui) {}

void godel_simple_gui::ExecutingState::onReset(BlendingWidget& gui) {}

void godel_simple_gui::ExecutingState::executeAll()
{
  try
  {
    for (std::size_t i = 0; i < plan_names_.size(); ++i)
    {
      execute_next_ = false;
      executeOne(plan_names_[i]);
      while(!execute_next_);
    }

    Q_EMIT newStateAvailable(new WaitToExecuteState(plan_names_));
  }
  catch (const BadExecutionError& err)
  {
    ROS_ERROR_STREAM("There was an error executing a plan " << err.what);
    Q_EMIT newStateAvailable(new ErrorState(err.what, new WaitToExecuteState(plan_names_)));
  }
}

void godel_simple_gui::ExecutingState::executeOne(const std::string& plan)
{
  ROS_DEBUG_STREAM("Executing " << plan);

  godel_msgs::SelectMotionPlanActionGoal goal;
  goal.goal.name = plan;
  goal.goal.simulate = false;
  goal.goal.wait_for_execution = true;
  select_motion_plan_action_client_.sendGoal(goal.goal,
    boost::bind(&godel_simple_gui::ExecutingState::selectMotionPlanDoneCallback, this, _1, _2));
}

void godel_simple_gui::ExecutingState::selectMotionPlanDoneCallback(const actionlib::SimpleClientGoalState& state,
            const godel_msgs::SelectMotionPlanResultConstPtr& result)
{
  execute_next_ = true;

  if (result->code != godel_msgs::SelectMotionPlanResult::SUCCESS)
  {
    std::ostringstream ss;
    ss << "Failed to execute plan. Please see logs for more details.";
    throw BadExecutionError(ss.str());
  }
}
