# Gemini Project Mandates

## Model Version Identification

**Mandate:** When the model version is required for a task (e.g., in a commit message), the following procedure must be followed to ensure accuracy:

1.  **State the need:** Clearly state that the model version is required for the upcoming action.
2.  **Ask the user:** Explicitly ask the user to provide the specific model version they are currently using.
3.  **Use the user-provided version:** Use the exact string provided by the user.

**Rationale:** The Gemini CLI agent does not have a built-in tool to programmatically access its own model version metadata from the underlying API. This mandate prevents the agent from making incorrect assumptions about its version and ensures that all reported model information is accurate and user-verified.
