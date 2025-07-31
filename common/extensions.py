"""
.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at
       http://www.apache.org/licenses/LICENSE-2.0
   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

import time

from ariadne.types import Extension, ContextValue


class QueryExecutionTimeExtension(Extension):
    def __init__(self):
        """Initialises the QueryExecutionTimeExtension instance.

        Sets the initial start timestamp to None.
        """
        self.start_timestamp = None

    def request_started(self, context: ContextValue) -> None:
        """Records the start timestamp when a request is initiated.

        Args:
            context (ContextValue): The context of the current request.
        """
        self.start_timestamp = time.perf_counter_ns()

    def format(self, context):
        """Formats and returns the execution time of the query.

        Calculates the execution time in seconds from the recorded start time.

        Args:
            context: The context of the current request.

        Returns:
            dict: A dictionary with the key 'execution_time_in_seconds' and its corresponding value.
        """
        if self.start_timestamp:
            exec_time_in_secs = round(
                (time.perf_counter_ns() - self.start_timestamp) / 1000000000, 2
            )
            return {"execution_time_in_seconds": exec_time_in_secs}
